from __future__ import annotations

import gzip

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pymatgen.core.tensors import Tensor
from pymatgen.io.aims.parsers import (
    EV_PER_A3_TO_KBAR,
    LINE_NOT_FOUND,
    AimsOutCalcChunk,
    AimsOutChunk,
    AimsOutHeaderChunk,
    AimsParseError,
)
from pymatgen.util.testing import TEST_FILES_DIR

PARSER_FILE_DIR = TEST_FILES_DIR / "io/aims/parser_checks"


EPS_HP = 1e-15  # The epsilon value used to compare numbers that are high-precision
EPS_LP = 1e-7  # The epsilon value used to compare numbers that are low-precision


@pytest.fixture
def default_chunk():
    lines = ["TEST", "A", "TEST", "| Number of atoms: 200 atoms"]
    return AimsOutChunk(lines)


def test_reverse_search_for(default_chunk):
    assert default_chunk.reverse_search_for(["TEST"]) == 2
    assert default_chunk.reverse_search_for(["TEST"], 1) == 2

    assert default_chunk.reverse_search_for(["TEST A"]) == LINE_NOT_FOUND

    assert default_chunk.reverse_search_for(["A"]) == 1
    assert default_chunk.reverse_search_for(["A"], 2) == LINE_NOT_FOUND


def test_search_for_all(default_chunk):
    assert default_chunk.search_for_all("TEST") == [0, 2]
    assert default_chunk.search_for_all("TEST", 0, 1) == [0]
    assert default_chunk.search_for_all("TEST", 1, -1) == [2]


def test_search_parse_scalar(default_chunk):
    assert default_chunk.parse_scalar("n_atoms") == 200
    assert default_chunk.parse_scalar("n_electrons") is None


@pytest.fixture
def empty_header_chunk():
    return AimsOutHeaderChunk([])


@pytest.mark.parametrize("attr_name", ["n_atoms", "n_bands", "n_electrons", "n_spins", "initial_structure"])
def test_missing_parameter(attr_name, empty_header_chunk):
    with pytest.raises(AimsParseError, match="No information about"):
        getattr(empty_header_chunk, attr_name)


def test_default_header_electronic_temperature(empty_header_chunk):
    assert empty_header_chunk.electronic_temperature == 0.0


def test_default_header_initial_lattice(empty_header_chunk):
    assert empty_header_chunk.initial_lattice is None


def test_default_header_is_md(empty_header_chunk):
    assert not empty_header_chunk.is_md


def test_default_header_is_relaxation(empty_header_chunk):
    assert not empty_header_chunk.is_relaxation


def test_default_header_n_k_points(empty_header_chunk):
    assert empty_header_chunk.n_k_points is None


def test_default_header_k_points(empty_header_chunk):
    assert empty_header_chunk.k_points is None


def test_default_header_k_point_weights(empty_header_chunk):
    assert empty_header_chunk.k_point_weights is None


@pytest.fixture
def initial_lattice():
    return np.array(
        [
            [1, 2.70300000, 3.70300000],
            [4.70300000, 2, 6.70300000],
            [8.70300000, 7.70300000, 3],
        ]
    )


@pytest.fixture
def header_chunk():
    with gzip.open(f"{PARSER_FILE_DIR}/header_chunk.out.gz", mode="rt") as hc_file:
        lines = hc_file.readlines()

    for ll, line in enumerate(lines):
        lines[ll] = line.strip()

    return AimsOutHeaderChunk(lines)


def test_header_n_atoms(header_chunk):
    assert header_chunk.n_atoms == 2


def test_header_n_bands(header_chunk):
    assert header_chunk.n_bands == 3


def test_header_n_electrons(header_chunk):
    assert header_chunk.n_electrons == 28


def test_header_n_spins(header_chunk):
    assert header_chunk.n_spins == 2


def test_header_initial_structure(header_chunk, initial_lattice):
    initial_positions = np.array([[0.000, 1.000, 2.000], [2.703, 3.703, 4.703]])
    assert len(header_chunk.initial_structure) == 2
    assert_allclose(header_chunk.initial_structure.lattice.matrix, initial_lattice)
    assert_allclose(header_chunk.initial_structure.cart_coords, initial_positions, atol=1e-12)
    assert [sp.symbol for sp in header_chunk.initial_structure.species] == ["Na", "Cl"]


def test_header_initial_lattice(header_chunk, initial_lattice):
    assert_allclose(header_chunk.initial_lattice.matrix, initial_lattice)


def test_header_electronic_temperature(header_chunk):
    assert header_chunk.electronic_temperature == 0.05


def test_header_is_md(header_chunk):
    assert header_chunk.is_md


def test_header_is_relaxation(header_chunk):
    assert header_chunk.is_relaxation


def test_header_n_k_points(header_chunk):
    assert header_chunk.n_k_points == 8


@pytest.fixture
def k_points():
    return np.array(
        [
            [0.000, 0.000, 0.000],
            [0.000, 0.000, 0.500],
            [0.000, 0.500, 0.000],
            [0.000, 0.500, 0.500],
            [0.500, 0.000, 0.000],
            [0.500, 0.000, 0.500],
            [0.500, 0.500, 0.000],
            [0.500, 0.500, 0.500],
        ]
    )


def test_header_k_point_weights(
    header_chunk,
):
    assert_allclose(header_chunk.k_point_weights, np.full((8), 0.125))


def test_header_k_points(header_chunk, k_points):
    assert_allclose(header_chunk.k_points, k_points)


def test_header_header_summary(header_chunk, k_points):
    header_summary = {
        "initial_structure": header_chunk.initial_structure,
        "initial_lattice": header_chunk.initial_lattice,
        "is_relaxation": True,
        "is_md": True,
        "n_atoms": 2,
        "n_bands": 3,
        "n_electrons": 28,
        "n_spins": 2,
        "electronic_temperature": 0.05,
        "n_k_points": 8,
        "k_points": k_points,
        "k_point_weights": np.full((8), 0.125),
    }
    for key, val in header_chunk.header_summary.items():
        if isinstance(val, np.ndarray):
            assert_allclose(val, header_summary[key])
        else:
            assert val == header_summary[key]


@pytest.fixture
def empty_calc_chunk(header_chunk):
    return AimsOutCalcChunk([], header_chunk)


def test_header_transfer_n_atoms(empty_calc_chunk):
    assert empty_calc_chunk.n_atoms == 2


def test_header_transfer_n_bands(empty_calc_chunk):
    assert empty_calc_chunk.n_bands == 3


def test_header_transfer_n_electrons(empty_calc_chunk):
    assert empty_calc_chunk.n_electrons == 28


def test_header_transfer_n_spins(empty_calc_chunk):
    assert empty_calc_chunk.n_spins == 2


def test_header_transfer_initial_lattice(empty_calc_chunk, initial_lattice):
    assert_allclose(empty_calc_chunk.initial_lattice.matrix, initial_lattice)


def test_header_transfer_initial_structure(empty_calc_chunk, initial_lattice):
    initial_positions = np.array([[0.000, 1.000, 2.000], [2.703, 3.703, 4.703]])

    assert_allclose(empty_calc_chunk.initial_structure.lattice.matrix, empty_calc_chunk.initial_lattice.matrix)
    assert len(empty_calc_chunk.initial_structure) == 2
    assert_allclose(empty_calc_chunk.initial_structure.lattice.matrix, initial_lattice)
    assert_allclose(empty_calc_chunk.initial_structure.cart_coords, initial_positions, atol=1e-12)
    assert [sp.symbol for sp in empty_calc_chunk.initial_structure.species] == ["Na", "Cl"]


def test_header_transfer_electronic_temperature(empty_calc_chunk):
    assert empty_calc_chunk.electronic_temperature == 0.05


def test_header_transfer_n_k_points(empty_calc_chunk):
    assert empty_calc_chunk.n_k_points == 8


def test_header_transfer_k_point_weights(empty_calc_chunk):
    assert_allclose(empty_calc_chunk.k_point_weights, np.full((8), 0.125))


def test_header_transfer_k_points(empty_calc_chunk, k_points):
    assert_allclose(empty_calc_chunk.k_points, k_points)


def test_default_calc_energy_raises_error(empty_calc_chunk):
    with pytest.raises(AimsParseError, match="No energy is associated with the structure."):
        _ = empty_calc_chunk.energy


@pytest.mark.parametrize(
    "attr",
    [
        "forces",
        "stresses",
        "stress",
        "free_energy",
        "n_iter",
        "magmom",
        "E_f",
        "dipole",
        "mulliken_charges",
        "mulliken_spins",
        "hirshfeld_charges",
        "hirshfeld_volumes",
        "hirshfeld_atomic_dipoles",
        "hirshfeld_dipole",
    ],
)
def test_chunk_defaults_none(attr, empty_calc_chunk):
    assert getattr(empty_calc_chunk, attr) is None


def test_default_calc_is_metallic(empty_calc_chunk):
    assert not empty_calc_chunk.is_metallic


def test_default_calc_converged(empty_calc_chunk):
    assert not empty_calc_chunk.converged


@pytest.fixture
def calc_chunk(header_chunk):
    with gzip.open(f"{PARSER_FILE_DIR}/calc_chunk.out.gz", mode="rt") as file:
        lines = file.readlines()

    for ll, line in enumerate(lines):
        lines[ll] = line.strip()
    return AimsOutCalcChunk(lines, header_chunk)


@pytest.fixture
def numerical_stress_chunk(header_chunk):
    with gzip.open(f"{PARSER_FILE_DIR}/numerical_stress.out.gz", mode="rt") as file:
        lines = file.readlines()

    for ll, line in enumerate(lines):
        lines[ll] = line.strip()
    return AimsOutCalcChunk(lines, header_chunk)


def test_calc_structure(calc_chunk, initial_lattice):
    initial_positions = np.array([[0.000, 1.000, 2.000], [2.703, 3.703, 4.703]])

    assert len(calc_chunk.structure.species) == 2
    assert_allclose(calc_chunk.structure.lattice.matrix, initial_lattice)
    assert_allclose(calc_chunk.structure.cart_coords, initial_positions, atol=1e-12)
    assert [sp.symbol for sp in calc_chunk.structure.species] == ["Na", "Cl"]


def test_calc_forces(calc_chunk):
    forces = np.array([[1.0, 2.0, 3.0], [6.0, 5.0, 4.0]])
    assert_allclose(calc_chunk.forces, forces)

    # Different because of the constraints
    assert_allclose(calc_chunk.structure.site_properties["force"], forces)
    assert_allclose(calc_chunk.results["forces"], forces)


def test_calc_stresses(calc_chunk):
    stresses = EV_PER_A3_TO_KBAR * np.array(
        [
            Tensor.from_voigt([-10.0, -20.0, -30.0, -60.0, -50.0, -40.0]),
            Tensor.from_voigt([10.0, 20.0, 30.0, 60.0, 50.0, 40.0]),
        ]
    )
    assert_allclose(calc_chunk.stresses, stresses)
    assert_allclose(calc_chunk.structure.site_properties["atomic_virial_stress"], stresses)
    assert_allclose(calc_chunk.results["stresses"], stresses)


def test_calc_stress(calc_chunk):
    stress = EV_PER_A3_TO_KBAR * np.array([[1.0, 2.0, 3.0], [2.0, 5.0, 6.0], [3.0, 6.0, 7.0]])
    assert_allclose(calc_chunk.stress, stress)
    assert_allclose(calc_chunk.structure.properties["stress"], stress)
    assert_allclose(calc_chunk.results["stress"], stress)


def test_calc_num_stress(numerical_stress_chunk):
    stress = EV_PER_A3_TO_KBAR * np.array([[1.0, 2.0, 3.0], [2.0, 5.0, 6.0], [3.0, 6.0, 7.0]])
    assert_allclose(numerical_stress_chunk.stress, stress)
    assert_allclose(numerical_stress_chunk.structure.properties["stress"], stress)
    assert_allclose(numerical_stress_chunk.results["stress"], stress)


def test_calc_free_energy(calc_chunk):
    free_energy = -3.169503986610555e05
    assert np.abs(calc_chunk.free_energy - free_energy) < EPS_HP
    assert np.abs(calc_chunk.structure.properties["free_energy"] - free_energy) < EPS_HP
    assert np.abs(calc_chunk.results["free_energy"] - free_energy) < EPS_HP


def test_calc_energy(calc_chunk):
    energy = -2.169503986610555e05
    assert np.abs(calc_chunk.energy - energy) < EPS_HP
    assert np.abs(calc_chunk.structure.properties["energy"] - energy) < EPS_HP
    assert np.abs(calc_chunk.results["energy"] - energy) < EPS_HP


def test_calc_magnetic_moment(calc_chunk):
    magmom = 0
    assert calc_chunk.magmom == magmom
    assert calc_chunk.structure.properties["magmom"] == magmom
    assert calc_chunk.results["magmom"] == magmom


def test_calc_n_iter(calc_chunk):
    n_iter = 58
    assert calc_chunk.n_iter == n_iter
    assert calc_chunk.results["n_iter"] == n_iter


def test_calc_fermi_energy(calc_chunk):
    Ef = -8.24271207
    assert np.abs(calc_chunk.E_f - Ef) < EPS_LP
    assert np.abs(calc_chunk.results["fermi_energy"] - Ef) < EPS_LP


def test_calc_dipole(calc_chunk):
    assert calc_chunk.dipole is None


def test_calc_is_metallic(calc_chunk):
    assert calc_chunk.is_metallic


def test_calc_converged(calc_chunk):
    assert calc_chunk.converged


def test_calc_mulliken_charges(calc_chunk):
    mulliken_charges = [0.617623, -0.617623]
    assert_allclose(calc_chunk.mulliken_charges, mulliken_charges)
    assert_allclose(calc_chunk.results["mulliken_charges"], mulliken_charges)


def test_calc_mulliken_spins(calc_chunk):
    # TARP: False numbers added to test parsing
    mulliken_spins = [-0.003141, 0.002718]
    assert_allclose(calc_chunk.mulliken_spins, mulliken_spins)
    assert_allclose(calc_chunk.results["mulliken_spins"], mulliken_spins)


def test_calc_hirshfeld_charges(calc_chunk):
    hirshfeld_charges = [0.20898543, -0.20840994]
    assert_allclose(calc_chunk.hirshfeld_charges, hirshfeld_charges)
    assert_allclose(calc_chunk.results["hirshfeld_charges"], hirshfeld_charges)


def test_calc_hirshfeld_volumes(calc_chunk):
    hirshfeld_volumes = [73.39467444, 62.86011074]
    assert_allclose(calc_chunk.hirshfeld_volumes, hirshfeld_volumes)
    assert_allclose(calc_chunk.results["hirshfeld_volumes"], hirshfeld_volumes)


def test_calc_hirshfeld_atomic_dipoles(calc_chunk):
    hirshfeld_atomic_dipoles = np.zeros((2, 3))
    assert_allclose(calc_chunk.hirshfeld_atomic_dipoles, hirshfeld_atomic_dipoles)
    assert_allclose(calc_chunk.results["hirshfeld_atomic_dipoles"], hirshfeld_atomic_dipoles)


def test_calc_hirshfeld_dipole(calc_chunk):
    assert calc_chunk.hirshfeld_dipole is None


@pytest.fixture
def molecular_header_chunk():
    with gzip.open(f"{PARSER_FILE_DIR}/molecular_header_chunk.out.gz", mode="rt") as file:
        lines = file.readlines()

    for ll, line in enumerate(lines):
        lines[ll] = line.strip()

    return AimsOutHeaderChunk(lines)


@pytest.mark.parametrize(
    "attr_name",
    ["k_points", "k_point_weights", "initial_lattice", "n_k_points"],
)
def test_chunk_molecular_header_defaults_none(attr_name, molecular_header_chunk):
    assert getattr(molecular_header_chunk, attr_name) is None


def test_molecular_header_n_bands(molecular_header_chunk):
    assert molecular_header_chunk.n_bands == 7


def test_molecular_header_initial_structure(molecular_header_chunk, molecular_positions):
    assert len(molecular_header_chunk.initial_structure) == 3
    assert [sp.symbol for sp in molecular_header_chunk.initial_structure.species] == ["O", "H", "H"]
    assert_allclose(
        molecular_header_chunk.initial_structure.cart_coords,
        [[0, 0, 0], [0.95840000, 0, 0], [-0.24000000, 0.92790000, 0]],
    )


@pytest.fixture
def molecular_calc_chunk(molecular_header_chunk):
    with gzip.open(f"{PARSER_FILE_DIR}/molecular_calc_chunk.out.gz", mode="rt") as file:
        lines = file.readlines()

    for idx, line in enumerate(lines):
        lines[idx] = line.strip()
    return AimsOutCalcChunk(lines, molecular_header_chunk)


@pytest.fixture
def molecular_positions():
    return np.array([[-0.00191785, -0.00243279, 0], [0.97071531, -0.00756333, 0], [-0.25039746, 0.93789612, 0]])


def test_molecular_calc_atoms(molecular_calc_chunk, molecular_positions):
    assert len(molecular_calc_chunk.structure.species) == 3
    assert_allclose(molecular_calc_chunk.structure.cart_coords, molecular_positions)
    assert [sp.symbol for sp in molecular_calc_chunk.structure.species] == ["O", "H", "H"]


def test_molecular_calc_forces(molecular_calc_chunk):
    forces = np.array(
        [
            [0.502371357164392e-03, 0.518627676606471e-03, 0.000000000000000e00],
            [-0.108826758257187e-03, -0.408128912334209e-03, -0.649037698626122e-27],
            [-0.393544598907207e-03, -0.110498764272267e-03, -0.973556547939183e-27],
        ]
    )
    assert_allclose(molecular_calc_chunk.forces, forces)
    assert_allclose(molecular_calc_chunk.structure.site_properties["force"], forces)
    assert_allclose(molecular_calc_chunk.results["forces"], forces)


@pytest.mark.parametrize("attrname", ["stresses", "stress", "magmom", "E_f"])
def test_chunk_molecular_defaults_none(attrname, molecular_calc_chunk):
    assert getattr(molecular_calc_chunk, attrname) is None


def test_molecular_calc_free_energy(molecular_calc_chunk):
    free_energy = -2.206778551123339e04
    assert np.abs(molecular_calc_chunk.free_energy - free_energy) < EPS_HP
    assert np.abs(molecular_calc_chunk.results["free_energy"] - free_energy) < EPS_HP
    assert np.abs(molecular_calc_chunk.structure.properties["free_energy"] - free_energy) < EPS_HP


def test_molecular_calc_energy(molecular_calc_chunk):
    energy = -0.206778551123339e04
    assert np.abs(molecular_calc_chunk.energy - energy) < EPS_HP
    assert np.abs(molecular_calc_chunk.structure.properties["energy"] - energy) < EPS_HP
    assert np.abs(molecular_calc_chunk.results["energy"] - energy) < EPS_HP


def test_molecular_calc_n_iter(molecular_calc_chunk):
    n_iter = 7
    assert molecular_calc_chunk.n_iter == n_iter
    assert molecular_calc_chunk.results["n_iter"] == n_iter


def test_molecular_calc_dipole(molecular_calc_chunk):
    dipole = [0.260286493869765, 0.336152447755231, 0.470003778119121e-15]
    assert_allclose(molecular_calc_chunk.dipole, dipole)
    assert_allclose(molecular_calc_chunk.structure.properties["dipole"], dipole)
    assert_allclose(molecular_calc_chunk.results["dipole"], dipole)


def test_molecular_calc_is_metallic(molecular_calc_chunk):
    assert not molecular_calc_chunk.is_metallic


def test_molecular_calc_converged(molecular_calc_chunk):
    assert molecular_calc_chunk.converged


def test_molecular_calc_hirshfeld_charges(molecular_calc_chunk):
    molecular_hirshfeld_charges = np.array([-0.32053200, 0.16022630, 0.16020375])
    assert_allclose(molecular_calc_chunk.hirshfeld_charges, molecular_hirshfeld_charges)
    assert_allclose(molecular_calc_chunk.results["hirshfeld_charges"], molecular_hirshfeld_charges)


def test_molecular_calc_hirshfeld_volumes(molecular_calc_chunk):
    hirshfeld_volumes = np.array([21.83060659, 6.07674041, 6.07684447])
    assert_allclose(molecular_calc_chunk.hirshfeld_volumes, hirshfeld_volumes)
    assert_allclose(molecular_calc_chunk.results["hirshfeld_volumes"], hirshfeld_volumes)


def test_molecular_calc_hirshfeld_atomic_dipoles(molecular_calc_chunk):
    hirshfeld_atomic_dipoles = np.array(
        [[0.04249319, 0.05486053, 0], [0.13710134, -0.00105126, 0], [-0.03534982, 0.13248706, 0]]
    )
    assert_allclose(molecular_calc_chunk.hirshfeld_atomic_dipoles, hirshfeld_atomic_dipoles)
    assert_allclose(
        molecular_calc_chunk.results["hirshfeld_atomic_dipoles"],
        hirshfeld_atomic_dipoles,
    )


def test_molecular_calc_hirshfeld_dipole(molecular_calc_chunk, molecular_positions):
    molecular_hirshfeld_charges = np.array([-0.32053200, 0.16022630, 0.16020375])
    hirshfeld_dipole = np.sum(molecular_hirshfeld_charges.reshape((-1, 1)) * molecular_positions, axis=1)

    assert_allclose(molecular_calc_chunk.hirshfeld_dipole, hirshfeld_dipole)
    assert_allclose(molecular_calc_chunk.results["hirshfeld_dipole"], hirshfeld_dipole)
