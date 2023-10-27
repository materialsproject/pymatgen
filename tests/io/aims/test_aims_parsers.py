# flake8: noqa
import numpy as np
from ase.io import read
from pymatgen.io.aims.parsers import (
    AimsParseError,
    AimsOutChunk,
    AimsOutHeaderChunk,
    AimsOutCalcChunk,
    EV_PER_A3_TO_KBAR,
    LINE_NOT_FOUND,
)

from ase.stress import full_3x3_to_voigt_6_stress

from numpy.linalg import norm

import pytest

eps_hp = 1e-15  # The espsilon value used to compare numbers that are high-precision
eps_lp = 1e-7  # The espsilon value used to compare numbers that are low-precision


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


@pytest.mark.parametrize(
    "attrname",
    [
        "n_atoms",
        "n_bands",
        "n_electrons",
        "n_spins",
        "initial_structure",
    ],
)
def test_missing_parameter(attrname, empty_header_chunk):
    with pytest.raises(AimsParseError, match="No information about"):
        getattr(empty_header_chunk, attrname)


def test_default_header_electronic_temperature(empty_header_chunk):
    assert empty_header_chunk.electronic_temperature == 0.1


# def test_default_header_constraints(empty_header_chunk):
#     assert empty_header_chunk.constraints == []


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
            [1.00000000, 2.70300000, 3.70300000],
            [4.70300000, 2.00000000, 6.70300000],
            [8.70300000, 7.70300000, 3.00000000],
        ]
    )


@pytest.fixture
def initial_positions():
    return (np.array([[0.000, 1.000, 2.000], [2.703, 3.703, 4.703]]),)


@pytest.fixture
def header_chunk():
    lines = """
        | Number of atoms                   :        2
        | Number of Kohn-Sham states (occupied + empty):       3
        The structure contains        2 atoms,  and a total of         28.000 electrons.
        | Number of spin channels           :        2
        Occupation type: Gaussian broadening, width =   0.500000E-01 eV.
        Found relaxation constraint for atom      1: x coordinate fixed.
        Found relaxation constraint for atom      1: y coordinate fixed.
        Found relaxation constraint for atom      2: All coordinates fixed.
        Input geometry:
        | Unit cell:
        |        1.00000000        2.70300000        3.70300000
        |        4.70300000        2.00000000        6.70300000
        |        8.70300000        7.70300000        3.00000000
        | Atomic structure:
        |       Atom                x [A]            y [A]            z [A]
        |    1: Species Na            0.00000000        1.00000000        2.00000000
        |    2: Species Cl            2.70300000        3.70300000        4.70300000
        Initializing the k-points
        Using symmetry for reducing the k-points
        | k-points reduced from:        8 to        8
        | Number of k-points                             :         8
        The eigenvectors in the calculations are REAL.
        | K-points in task   0:         2
        | K-points in task   1:         2
        | K-points in task   2:         2
        | K-points in task   3:         2
        | k-point: 1 at     0.000000    0.000000    0.000000  , weight:   0.12500000
        | k-point: 2 at     0.000000    0.000000    0.500000  , weight:   0.12500000
        | k-point: 3 at     0.000000    0.500000    0.000000  , weight:   0.12500000
        | k-point: 4 at     0.000000    0.500000    0.500000  , weight:   0.12500000
        | k-point: 5 at     0.500000    0.000000    0.000000  , weight:   0.12500000
        | k-point: 6 at     0.500000    0.000000    0.500000  , weight:   0.12500000
        | k-point: 7 at     0.500000    0.500000    0.000000  , weight:   0.12500000
        | k-point: 8 at     0.500000    0.500000    0.500000  , weight:   0.12500000
        Geometry relaxation: A file "geometry.in.next_step" is written out by default after each step.
        Complete information for previous time-step:
    """
    lines = lines.splitlines()
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


def test_header_initial_structure(header_chunk, initial_lattice, initial_positions):
    assert len(header_chunk.initial_structure) == 2
    assert np.allclose(
        header_chunk.initial_structure.lattice.matrix,
        initial_lattice,
    )
    assert np.allclose(header_chunk.initial_structure.cart_coords, initial_positions)
    assert np.all(["Na", "Cl"] == [sp.symbol for sp in header_chunk.initial_structure.species])


def test_header_initial_lattice(header_chunk, initial_lattice):
    assert np.allclose(header_chunk.initial_lattice.matrix, initial_lattice)


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


@pytest.fixture
def k_point_weights():
    return np.full((8), 0.125)


def test_header_k_point_weights(header_chunk, k_point_weights):
    assert np.allclose(header_chunk.k_point_weights, k_point_weights)


def test_header_k_points(header_chunk, k_points):
    assert np.allclose(header_chunk.k_points, k_points)


def test_header_header_summary(header_chunk, k_points, k_point_weights):
    header_summary = {
        "initial_structure": header_chunk.initial_structure,
        "initial_lattice": header_chunk.initial_lattice,
        # "constraints": header_chunk.constraints,
        "is_relaxation": True,
        "is_md": True,
        "n_atoms": 2,
        "n_bands": 3,
        "n_electrons": 28,
        "n_spins": 2,
        "electronic_temperature": 0.05,
        "n_k_points": 8,
        "k_points": k_points,
        "k_point_weights": k_point_weights,
    }
    for key, val in header_chunk.header_summary.items():
        # if key == "constraints":
        #     continue
        if isinstance(val, np.ndarray):
            assert np.allclose(val, header_summary[key])
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
    assert np.allclose(empty_calc_chunk.initial_lattice.matrix, initial_lattice)


def test_header_transfer_initial_structure(empty_calc_chunk, initial_lattice, initial_positions):
    assert np.allclose(empty_calc_chunk.initial_structure.lattice.matrix, empty_calc_chunk.initial_lattice.matrix)
    assert len(empty_calc_chunk.initial_structure) == 2
    assert np.allclose(empty_calc_chunk.initial_structure.lattice.matrix, initial_lattice)
    assert np.allclose(empty_calc_chunk.initial_structure.cart_coords, initial_positions)
    assert np.all(["Na", "Cl"] == [sp.symbol for sp in empty_calc_chunk.initial_structure.species])


def test_header_transfer_electronic_temperature(empty_calc_chunk):
    assert empty_calc_chunk.electronic_temperature == 0.05


def test_header_transfer_n_k_points(empty_calc_chunk):
    assert empty_calc_chunk.n_k_points == 8


def test_header_transfer_k_point_weights(empty_calc_chunk, k_point_weights):
    assert np.allclose(empty_calc_chunk.k_point_weights, k_point_weights)


def test_header_transfer_k_points(empty_calc_chunk, k_points):
    assert np.allclose(empty_calc_chunk.k_points, k_points)


def test_default_calc_energy_raises_error(empty_calc_chunk):
    with pytest.raises(AimsParseError, match="No energy is associated with the structure."):
        getattr(empty_calc_chunk, "energy")


@pytest.mark.parametrize(
    "attrname",
    [
        "forces",
        "stresses",
        "stress",
        "free_energy",
        "n_iter",
        "magmom",
        "E_f",
        "dipole",
        "hirshfeld_charges",
        "hirshfeld_volumes",
        "hirshfeld_atomic_dipoles",
        "hirshfeld_dipole",
    ],
)
def test_chunk_defaults_none(attrname, empty_calc_chunk):
    assert getattr(empty_calc_chunk, attrname) is None


def test_default_calc_is_metallic(empty_calc_chunk):
    assert not empty_calc_chunk.is_metallic


def test_default_calc_converged(empty_calc_chunk):
    assert not empty_calc_chunk.converged


@pytest.fixture
def calc_chunk(header_chunk):
    lines = """
        | Number of self-consistency cycles          :           58
        | N = N_up - N_down (sum over all k points):         0.00000
        | Chemical potential (Fermi level):    -8.24271207 eV
        Total atomic forces (unitary forces were cleaned, then relaxation constraints were applied) [eV/Ang]:
        |    1          1.000000000000000E+00          2.000000000000000E+00          3.000000000000000E+00
        |    2          6.000000000000000E+00          5.000000000000000E+00          4.000000000000000E+00
        - Per atom stress (eV) used for heat flux calculation:
        Atom   | Stress components (1,1), (2,2), (3,3), (1,2), (1,3), (2,3)
        -------------------------------------------------------------------
        1 |    -1.0000000000E+01   -2.0000000000E+01   -3.0000000000E+01   -4.0000000000E+01   -5.0000000000E+01   -6.0000000000E+01
        2 |     1.0000000000E+01    2.0000000000E+01    3.0000000000E+01    4.0000000000E+01    5.0000000000E+01    6.0000000000E+01
        -------------------------------------------------------------------
        +-------------------------------------------------------------------+
        |              Analytical stress tensor - Symmetrized               |
        |                  Cartesian components [eV/A**3]                   |
        +-------------------------------------------------------------------+
        |                x                y                z                |
        |                                                                   |
        |  x         1.00000000       2.00000000       3.00000000           |
        |  y         2.00000000       5.00000000       6.00000000           |
        |  z         3.00000000       6.00000000       7.00000000           |
        |                                                                   |
        |  Pressure:       0.00383825   [eV/A**3]                           |
        |                                                                   |
        +-------------------------------------------------------------------+
        Energy and forces in a compact form:
        | Total energy uncorrected      :         -0.169503986610555E+05 eV
        | Total energy corrected        :         -2.169503986610555E+05 eV  <-- do not rely on this value for anything but (periodic) metals
        | Electronic free energy        :         -3.169503986610555E+05 eV
        Performing Hirshfeld analysis of fragment charges and moments.
        ----------------------------------------------------------------------
        | Atom     1: Na
        |   Hirshfeld charge        :      0.20898543
        |   Free atom volume        :     82.26734086
        |   Hirshfeld volume        :     73.39467444
        |   Hirshfeld spin moment   :      0.00000000
        |   Hirshfeld dipole vector :      0.00000000       0.00000000      -0.00000000
        |   Hirshfeld dipole moment :      0.00000000
        |   Hirshfeld second moments:      0.19955428       0.00000002       0.00000002
        |                                  0.00000002       0.19955473      -0.00000005
        |                                  0.00000002      -0.00000005       0.19955473
        ----------------------------------------------------------------------
        | Atom     2: Cl
        |   Hirshfeld charge        :     -0.20840994
        |   Free atom volume        :     65.62593663
        |   Hirshfeld volume        :     62.86011074
        |   Hirshfeld spin moment   :      0.00000000
        |   Hirshfeld dipole vector :      0.00000000       0.00000000       0.00000000
        |   Hirshfeld dipole moment :      0.00000000
        |   Hirshfeld second moments:      0.01482616      -0.00000001      -0.00000001
        |                                 -0.00000001       0.01482641       0.00000001
        |                                 -0.00000001       0.00000001       0.01482641
        ----------------------------------------------------------------------
        Writing Kohn-Sham eigenvalues.

        Spin-up eigenvalues:
        K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.036749           1.000000
        2        1             0.183747           5.000000
        3        0             0.330744           9.000000

        Spin-down eigenvalues:
        K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.110248           3.000000
        2        1             0.257245           7.000000
        3        0             0.404243           11.000000


        Spin-up eigenvalues:
        K-point:       2 at    0.000000    0.000000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.477741           13.000000
        2        1             0.624739           17.000000
        3        0             0.771736           21.000000

        Spin-down eigenvalues:
        K-point:       2 at    0.000000    0.000000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.551240           15.000000
        2        1             0.698237           19.000000
        3        0             0.845234           23.000000


        Spin-up eigenvalues:
        K-point:       3 at    0.000000    0.500000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.918733           25.000000
        2        1             1.065730           29.000000
        3        0             1.212728           33.000000

        Spin-down eigenvalues:
        K-point:       3 at    0.000000    0.500000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.992232           27.000000
        2        1             1.139229           31.000000
        3        0             1.286226           35.000000


        Spin-up eigenvalues:
        K-point:       4 at    0.000000    0.500000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             1.359725           37.000000
        2        1             1.506722           41.000000
        3        0             1.653720           45.000000

        Spin-down eigenvalues:
        K-point:       4 at    0.000000    0.500000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             1.433224           39.000000
        2        1             1.580221           43.000000
        3        0             1.727218           47.000000


        Spin-up eigenvalues:
        K-point:       5 at    0.500000    0.000000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             1.800717           49.000000
        2        1             1.947714           53.000000
        3        0             2.094711           57.000000

        Spin-down eigenvalues:
        K-point:       5 at    0.500000    0.000000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             1.874216           51.000000
        2        1             2.021213           55.000000
        3        0             2.168210           59.000000


        Spin-up eigenvalues:
        K-point:       6 at    0.500000    0.000000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             2.241709           61.000000
        2        1             2.388706           65.000000
        3        0             2.535703           69.000000

        Spin-down eigenvalues:
        K-point:       6 at    0.500000    0.000000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             2.315207           63.000000
        2        1             2.462205           67.000000
        3        0             2.609202           71.000000


        Spin-up eigenvalues:
        K-point:       7 at    0.500000    0.500000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             2.682701           73.000000
        2        1             2.829698           77.000000
        3        0             2.976695           81.000000

        Spin-down eigenvalues:
        K-point:       7 at    0.500000    0.500000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             2.756199           75.000000
        2        1             2.903197           79.000000
        3        0             3.050194           83.000000


        Spin-up eigenvalues:
        K-point:       8 at    0.500000    0.500000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             3.123693           85.000000
        2        1             3.270690           89.000000
        3        0             3.417687           93.000000

        Spin-down eigenvalues:
        K-point:       8 at    0.500000    0.500000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             3.197191           87.000000
        2        1             3.344189           91.000000
        3        0             3.491186           95.000000

        Current spin moment of the entire structure :
        | N = N_up - N_down (sum over all k points):         0.00000
        | S (sum over all k points)                :         0.00000

        What follows are estimated values for band gap, HOMO, LUMO, etc.
        | They are estimated on a discrete k-point grid and not necessarily exact.
        | For converged numbers, create a DOS and/or band structure plot on a denser k-grid.

        Highest occupied state (VBM) at     -8.19345940 eV (relative to internal zero)
        | Occupation number:      1.00000000
        | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)
        | Spin channel:        1

        Lowest unoccupied state (CBM) at    -3.62542909 eV (relative to internal zero)
        | Occupation number:      0.00000000
        | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)
        | Spin channel:        1

        ESTIMATED overall HOMO-LUMO gap:      4.56803031 eV between HOMO at k-point 1 and LUMO at k-point 1
        | This appears to be a direct band gap.
        The gap value is above 0.2 eV. Unless you are using a very sparse k-point grid
        this system is most likely an insulator or a semiconductor.
        | Chemical Potential                          :    -7.44914181 eV

        Self-consistency cycle converged.
        material is metallic within the approximate finite broadening function (occupation_type)
        Have a nice day.
        ------------------------------------------------------------
    """
    lines = lines.splitlines()
    for ll, line in enumerate(lines):
        lines[ll] = line.strip()
    return AimsOutCalcChunk(lines, header_chunk)


@pytest.fixture
def numerical_stress_chunk(header_chunk):
    lines = """
        | Number of self-consistency cycles          :           58
        | N = N_up - N_down (sum over all k points):         0.00000
        | Chemical potential (Fermi level):    -8.24271207 eV
        Total atomic forces (unitary forces were cleaned, then relaxation constraints were applied) [eV/Ang]:
        |    1          1.000000000000000E+00          2.000000000000000E+00          3.000000000000000E+00
        |    2          6.000000000000000E+00          5.000000000000000E+00          4.000000000000000E+00
        - Per atom stress (eV) used for heat flux calculation:
        Atom   | Stress components (1,1), (2,2), (3,3), (1,2), (1,3), (2,3)
        -------------------------------------------------------------------
        1 |    -1.0000000000E+01   -2.0000000000E+01   -3.0000000000E+01   -4.0000000000E+01   -5.0000000000E+01   -6.0000000000E+01
        2 |     1.0000000000E+01    2.0000000000E+01    3.0000000000E+01    4.0000000000E+01    5.0000000000E+01    6.0000000000E+01
        -------------------------------------------------------------------
        +-------------------------------------------------------------------+
        |                       Numerical stress tensor                     |
        |                    Cartesian components [eV/A**3]                 |
        +-------------------------------------------------------------------+
        |                x                y                z                |
        |                                                                   |
        |  x         1.00000000       2.00000000       3.00000000           |
        |  y         2.00000000       5.00000000       6.00000000           |
        |  z         3.00000000       6.00000000       7.00000000           |
        |                                                                   |
        |  Pressure:       0.00383825   [eV/A**3]                           |
        |                                                                   |
        +-------------------------------------------------------------------+
        Energy and forces in a compact form:
        | Total energy uncorrected      :         -0.169503986610555E+05 eV
        | Total energy corrected        :         -2.169503986610555E+05 eV  <-- do not rely on this value for anything but (periodic) metals
        | Electronic free energy        :         -3.169503986610555E+05 eV

        Writing Kohn-Sham eigenvalues.

        Spin-up eigenvalues:
        K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.036749           1.000000
        2        1             0.183747           5.000000
        3        0             0.330744           9.000000

        Spin-down eigenvalues:
        K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.110248           3.000000
        2        1             0.257245           7.000000
        3        0             0.404243           11.000000


        Spin-up eigenvalues:
        K-point:       2 at    0.000000    0.000000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.477741           13.000000
        2        1             0.624739           17.000000
        3        0             0.771736           21.000000

        Spin-down eigenvalues:
        K-point:       2 at    0.000000    0.000000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.551240           15.000000
        2        1             0.698237           19.000000
        3        0             0.845234           23.000000


        Spin-up eigenvalues:
        K-point:       3 at    0.000000    0.500000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.918733           25.000000
        2        1             1.065730           29.000000
        3        0             1.212728           33.000000

        Spin-down eigenvalues:
        K-point:       3 at    0.000000    0.500000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             0.992232           27.000000
        2        1             1.139229           31.000000
        3        0             1.286226           35.000000


        Spin-up eigenvalues:
        K-point:       4 at    0.000000    0.500000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             1.359725           37.000000
        2        1             1.506722           41.000000
        3        0             1.653720           45.000000

        Spin-down eigenvalues:
        K-point:       4 at    0.000000    0.500000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             1.433224           39.000000
        2        1             1.580221           43.000000
        3        0             1.727218           47.000000


        Spin-up eigenvalues:
        K-point:       5 at    0.500000    0.000000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             1.800717           49.000000
        2        1             1.947714           53.000000
        3        0             2.094711           57.000000

        Spin-down eigenvalues:
        K-point:       5 at    0.500000    0.000000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             1.874216           51.000000
        2        1             2.021213           55.000000
        3        0             2.168210           59.000000


        Spin-up eigenvalues:
        K-point:       6 at    0.500000    0.000000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             2.241709           61.000000
        2        1             2.388706           65.000000
        3        0             2.535703           69.000000

        Spin-down eigenvalues:
        K-point:       6 at    0.500000    0.000000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             2.315207           63.000000
        2        1             2.462205           67.000000
        3        0             2.609202           71.000000


        Spin-up eigenvalues:
        K-point:       7 at    0.500000    0.500000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             2.682701           73.000000
        2        1             2.829698           77.000000
        3        0             2.976695           81.000000

        Spin-down eigenvalues:
        K-point:       7 at    0.500000    0.500000    0.000000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             2.756199           75.000000
        2        1             2.903197           79.000000
        3        0             3.050194           83.000000


        Spin-up eigenvalues:
        K-point:       8 at    0.500000    0.500000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             3.123693           85.000000
        2        1             3.270690           89.000000
        3        0             3.417687           93.000000

        Spin-down eigenvalues:
        K-point:       8 at    0.500000    0.500000    0.500000 (in units of recip. lattice)

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1        1             3.197191           87.000000
        2        1             3.344189           91.000000
        3        0             3.491186           95.000000

        Current spin moment of the entire structure :
        | N = N_up - N_down (sum over all k points):         0.00000
        | S (sum over all k points)                :         0.00000

        What follows are estimated values for band gap, HOMO, LUMO, etc.
        | They are estimated on a discrete k-point grid and not necessarily exact.
        | For converged numbers, create a DOS and/or band structure plot on a denser k-grid.

        Highest occupied state (VBM) at     -8.19345940 eV (relative to internal zero)
        | Occupation number:      1.00000000
        | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)
        | Spin channel:        1

        Lowest unoccupied state (CBM) at    -3.62542909 eV (relative to internal zero)
        | Occupation number:      0.00000000
        | K-point:       1 at    0.000000    0.000000    0.000000 (in units of recip. lattice)
        | Spin channel:        1

        ESTIMATED overall HOMO-LUMO gap:      4.56803031 eV between HOMO at k-point 1 and LUMO at k-point 1
        | This appears to be a direct band gap.
        The gap value is above 0.2 eV. Unless you are using a very sparse k-point grid
        this system is most likely an insulator or a semiconductor.
        | Chemical Potential                          :    -7.44914181 eV

        Self-consistency cycle converged.
        material is metallic within the approximate finite broadening function (occupation_type)
        Have a nice day.
        ------------------------------------------------------------
    """
    lines = lines.splitlines()
    for ll, line in enumerate(lines):
        lines[ll] = line.strip()
    return AimsOutCalcChunk(lines, header_chunk)


def test_calc_structure(calc_chunk, initial_lattice, initial_positions):
    assert len(calc_chunk.structure.species) == 2
    assert np.allclose(calc_chunk.structure.lattice.matrix, initial_lattice)
    assert np.allclose(calc_chunk.structure.cart_coords, initial_positions)
    assert np.all(["Na", "Cl"] == [sp.symbol for sp in calc_chunk.structure.species])


def test_calc_forces(calc_chunk):
    forces = np.array([[1.0, 2.0, 3.0], [6.0, 5.0, 4.0]])
    assert np.allclose(calc_chunk.forces, forces)

    # Different because of the constraints
    assert np.allclose(calc_chunk.structure.site_properties["force"], forces)
    assert np.allclose(calc_chunk.results["forces"], forces)


def test_calc_stresses(calc_chunk):
    stresses = EV_PER_A3_TO_KBAR * np.array(
        [
            [-10.0, -20.0, -30.0, -60.0, -50.0, -40.0],
            [10.0, 20.0, 30.0, 60.0, 50.0, 40.0],
        ]
    )
    assert np.allclose(calc_chunk.stresses, stresses)
    assert np.allclose(calc_chunk.structure.site_properties["atomic_virial_stress"], stresses)
    assert np.allclose(calc_chunk.results["stresses"], stresses)


def test_calc_stress(calc_chunk):
    stress = EV_PER_A3_TO_KBAR * full_3x3_to_voigt_6_stress(
        np.array(
            [
                [1.00000000, 2.00000000, 3.00000000],
                [2.00000000, 5.00000000, 6.00000000],
                [3.00000000, 6.00000000, 7.00000000],
            ]
        )
    )
    assert np.allclose(calc_chunk.stress, stress)
    assert np.allclose(calc_chunk.structure.properties["stress"], stress)
    assert np.allclose(calc_chunk.results["stress"], stress)


def test_calc_num_stress(numerical_stress_chunk):
    stress = EV_PER_A3_TO_KBAR * full_3x3_to_voigt_6_stress(
        np.array(
            [
                [1.00000000, 2.00000000, 3.00000000],
                [2.00000000, 5.00000000, 6.00000000],
                [3.00000000, 6.00000000, 7.00000000],
            ]
        )
    )
    assert np.allclose(numerical_stress_chunk.stress, stress)
    assert np.allclose(numerical_stress_chunk.structure.properties["stress"], stress)
    assert np.allclose(numerical_stress_chunk.results["stress"], stress)


def test_calc_free_energy(calc_chunk):
    free_energy = -3.169503986610555e05
    assert np.abs(calc_chunk.free_energy - free_energy) < eps_hp
    assert np.abs(calc_chunk.structure.properties["free_energy"] - free_energy) < eps_hp
    assert np.abs(calc_chunk.results["free_energy"] - free_energy) < eps_hp


def test_calc_energy(calc_chunk):
    energy = -2.169503986610555e05
    assert np.abs(calc_chunk.energy - energy) < eps_hp
    assert np.abs(calc_chunk.structure.properties["energy"] - energy) < eps_hp
    assert np.abs(calc_chunk.results["energy"] - energy) < eps_hp


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
    assert np.abs(calc_chunk.E_f - Ef) < eps_lp
    assert np.abs(calc_chunk.results["fermi_energy"] - Ef) < eps_lp


def test_calc_dipole(calc_chunk):
    assert calc_chunk.dipole is None


def test_calc_is_metallic(calc_chunk):
    assert calc_chunk.is_metallic


def test_calc_converged(calc_chunk):
    assert calc_chunk.converged


def test_calc_hirshfeld_charges(calc_chunk):
    hirshfeld_charges = [0.20898543, -0.20840994]
    assert np.allclose(calc_chunk.hirshfeld_charges, hirshfeld_charges)
    assert np.allclose(calc_chunk.results["hirshfeld_charges"], hirshfeld_charges)


def test_calc_hirshfeld_volumes(calc_chunk):
    hirshfeld_volumes = [73.39467444, 62.86011074]
    assert np.allclose(calc_chunk.hirshfeld_volumes, hirshfeld_volumes)
    assert np.allclose(calc_chunk.results["hirshfeld_volumes"], hirshfeld_volumes)


def test_calc_hirshfeld_atomic_dipoles(calc_chunk):
    hirshfeld_atomic_dipoles = np.zeros((2, 3))
    assert np.allclose(calc_chunk.hirshfeld_atomic_dipoles, hirshfeld_atomic_dipoles)
    assert np.allclose(calc_chunk.results["hirshfeld_atomic_dipoles"], hirshfeld_atomic_dipoles)


def test_calc_hirshfeld_dipole(calc_chunk):
    assert calc_chunk.hirshfeld_dipole is None


@pytest.fixture
def molecular_header_chunk():
    lines = """
        | Number of atoms                   :        3
        | Number of spin channels           :        1
        The structure contains        3 atoms  and a total of         10.000 electrons.
        Input geometry:
        | Atomic structure:
        |       Atom                x [A]            y [A]            z [A]
        |    1: Species O             0.00000000        0.00000000        0.00000000
        |    2: Species H             0.95840000        0.00000000        0.00000000
        |    3: Species H            -0.24000000        0.92790000        0.00000000
        'Geometry relaxation: A file geometry.in.next_step is written out by default after each step.'
        | Maximum number of basis functions            :        7
        | Number of Kohn-Sham states (occupied + empty):       11
        Reducing total number of  Kohn-Sham states to        7.
    """

    lines = lines.splitlines()
    for ll, line in enumerate(lines):
        lines[ll] = line.strip()

    return AimsOutHeaderChunk(lines)


@pytest.mark.parametrize(
    "attrname",
    [
        "k_points",
        "k_point_weights",
        "initial_lattice",
        "n_k_points",
    ],
)
def test_chunk_molecular_header_defaults_none(attrname, molecular_header_chunk):
    assert getattr(molecular_header_chunk, attrname) is None


def test_molecular_header_n_bands(molecular_header_chunk):
    assert molecular_header_chunk.n_bands == 7


def test_molecular_header_initial_structure(molecular_header_chunk, molecular_positions):
    assert len(molecular_header_chunk.initial_structure) == 3
    assert np.all(["O", "H", "H"] == [sp.symbol for sp in molecular_header_chunk.initial_structure.species])
    assert np.allclose(
        molecular_header_chunk.initial_structure.cart_coords,
        np.array(
            [
                [0.00000000, 0.00000000, 0.00000000],
                [0.95840000, 0.00000000, 0.00000000],
                [-0.24000000, 0.92790000, 0.00000000],
            ]
        ),
    )


@pytest.fixture
def molecular_calc_chunk(molecular_header_chunk):
    lines = """
        | Number of self-consistency cycles          :           7
        | Chemical Potential                          :    -0.61315483 eV
        Updated atomic structure:
        x [A]             y [A]             z [A]
        atom        -0.00191785       -0.00243279        0.00000000  O
        atom         0.97071531       -0.00756333        0.00000000  H
        atom        -0.25039746        0.93789612       -0.00000000  H
        | Total dipole moment [eAng]          :          0.260286493869765E+00         0.336152447755231E+00         0.470003778119121E-15
        Energy and forces in a compact form:
        | Total energy uncorrected      :         -0.206778551123339E+04 eV
        | Total energy corrected        :         -5.206778551123339E+04 eV  <-- do not rely on this value for anything but (periodic) metals
        | Electronic free energy        :         -2.206778551123339E+04 eV
        Total atomic forces (unitary forces cleaned) [eV/Ang]:
        |    1          0.502371357164392E-03          0.518627676606471E-03          0.000000000000000E+00
        |    2         -0.108826758257187E-03         -0.408128912334209E-03         -0.649037698626122E-27
        |    3         -0.393544598907207E-03         -0.110498764272267E-03         -0.973556547939183E-27
        Performing Hirshfeld analysis of fragment charges and moments.
        ----------------------------------------------------------------------
        | Atom     1: O
        |   Hirshfeld charge        :     -0.32053200
        |   Free atom volume        :     23.59848617
        |   Hirshfeld volume        :     21.83060659
        |   Hirshfeld dipole vector :      0.04249319       0.05486053       0.00000000
        |   Hirshfeld dipole moment :      0.06939271
        |   Hirshfeld second moments:      0.04964380      -0.04453278      -0.00000000
        |                                 -0.04453278       0.02659295       0.00000000
        |                                 -0.00000000       0.00000000      -0.05608173
        ----------------------------------------------------------------------
        | Atom     2: H
        |   Hirshfeld charge        :      0.16022630
        |   Free atom volume        :     10.48483941
        |   Hirshfeld volume        :      6.07674041
        |   Hirshfeld dipole vector :      0.13710134      -0.00105126       0.00000000
        |   Hirshfeld dipole moment :      0.13710537
        |   Hirshfeld second moments:      0.12058896      -0.01198026      -0.00000000
        |                                 -0.01198026       0.14550360       0.00000000
        |                                 -0.00000000       0.00000000       0.10836357
        ----------------------------------------------------------------------
        | Atom     3: H
        |   Hirshfeld charge        :      0.16020375
        |   Free atom volume        :     10.48483941
        |   Hirshfeld volume        :      6.07684447
        |   Hirshfeld dipole vector :     -0.03534982       0.13248706       0.00000000
        |   Hirshfeld dipole moment :      0.13712195
        |   Hirshfeld second moments:      0.14974686      -0.00443579      -0.00000000
        |                                 -0.00443579       0.11633028      -0.00000000
        |                                 -0.00000000      -0.00000000       0.10836209
        ----------------
        Writing Kohn-Sham eigenvalues.

        State    Occupation    Eigenvalue [Ha]    Eigenvalue [eV]
        1       2.00000         -18.640915         -507.24511
        2       2.00000          -0.918449          -24.99226
        3       2.00000          -0.482216          -13.12175
        4       2.00000          -0.338691           -9.21626
        5       2.00000          -0.264427           -7.19543
        6       0.00000          -0.000414           -0.01127
        7       0.00000           0.095040            2.58616

        Highest occupied state (VBM) at     -7.19542820 eV
        | Occupation number:      2.00000000

        Lowest unoccupied state (CBM) at    -0.01126981 eV
        | Occupation number:      0.00000000

        Overall HOMO-LUMO gap:      7.18415839 eV.
        | Chemical Potential                          :    -0.61315483 eV

        Self-consistency cycle converged.
        Have a nice day.
        ------------------------------------------------------------

    """
    lines = lines.splitlines()
    for ll, line in enumerate(lines):
        lines[ll] = line.strip()
    return AimsOutCalcChunk(lines, molecular_header_chunk)


@pytest.fixture
def molecular_positions():
    return np.array(
        [
            [-0.00191785, -0.00243279, 0.00000000],
            [0.97071531, -0.00756333, 0.00000000],
            [-0.25039746, 0.93789612, 0.00000000],
        ]
    )


def test_molecular_calc_atoms(molecular_calc_chunk, molecular_positions):
    assert len(molecular_calc_chunk.structure.species) == 3
    assert np.allclose(molecular_calc_chunk.structure.cart_coords, molecular_positions)
    assert np.all(["O", "H", "H"] == [sp.symbol for sp in molecular_calc_chunk.structure.species])


def test_molecular_calc_forces(molecular_calc_chunk):
    forces = np.array(
        [
            [0.502371357164392e-03, 0.518627676606471e-03, 0.000000000000000e00],
            [-0.108826758257187e-03, -0.408128912334209e-03, -0.649037698626122e-27],
            [-0.393544598907207e-03, -0.110498764272267e-03, -0.973556547939183e-27],
        ]
    )
    assert np.allclose(molecular_calc_chunk.forces, forces)
    assert np.allclose(molecular_calc_chunk.structure.site_properties["force"], forces)
    assert np.allclose(molecular_calc_chunk.results["forces"], forces)


@pytest.mark.parametrize(
    "attrname",
    [
        "stresses",
        "stress",
        "magmom",
        "E_f",
    ],
)
def test_chunk_molecular_defaults_none(attrname, molecular_calc_chunk):
    assert getattr(molecular_calc_chunk, attrname) is None


def test_molecular_calc_free_energy(molecular_calc_chunk):
    free_energy = -2.206778551123339e04
    assert np.abs(molecular_calc_chunk.free_energy - free_energy) < eps_hp
    assert np.abs(molecular_calc_chunk.results["free_energy"] - free_energy) < eps_hp
    assert np.abs(molecular_calc_chunk.structure.properties["free_energy"] - free_energy) < eps_hp


def test_molecular_calc_energy(molecular_calc_chunk):
    energy = -0.206778551123339e04
    assert np.abs(molecular_calc_chunk.energy - energy) < eps_hp
    assert np.abs(molecular_calc_chunk.structure.properties["energy"] - energy) < eps_hp
    assert np.abs(molecular_calc_chunk.results["energy"] - energy) < eps_hp


def test_molecular_calc_n_iter(molecular_calc_chunk):
    n_iter = 7
    assert molecular_calc_chunk.n_iter == n_iter
    assert molecular_calc_chunk.results["n_iter"] == n_iter


def test_molecular_calc_dipole(molecular_calc_chunk):
    dipole = [0.260286493869765, 0.336152447755231, 0.470003778119121e-15]
    assert np.allclose(molecular_calc_chunk.dipole, dipole)
    assert np.allclose(molecular_calc_chunk.structure.properties["dipole"], dipole)
    assert np.allclose(molecular_calc_chunk.results["dipole"], dipole)


def test_molecular_calc_is_metallic(molecular_calc_chunk):
    assert not molecular_calc_chunk.is_metallic


def test_molecular_calc_converged(molecular_calc_chunk):
    assert molecular_calc_chunk.converged


@pytest.fixture
def molecular_hirshfeld_charges():
    return np.array([-0.32053200, 0.16022630, 0.16020375])


def test_molecular_calc_hirshfeld_charges(molecular_calc_chunk, molecular_hirshfeld_charges):
    assert np.allclose(molecular_calc_chunk.hirshfeld_charges, molecular_hirshfeld_charges)
    assert np.allclose(molecular_calc_chunk.results["hirshfeld_charges"], molecular_hirshfeld_charges)


def test_molecular_calc_hirshfeld_volumes(molecular_calc_chunk):
    hirshfeld_volumes = np.array([21.83060659, 6.07674041, 6.07684447])
    assert np.allclose(molecular_calc_chunk.hirshfeld_volumes, hirshfeld_volumes)
    assert np.allclose(molecular_calc_chunk.results["hirshfeld_volumes"], hirshfeld_volumes)


def test_molecular_calc_hirshfeld_atomic_dipoles(molecular_calc_chunk):
    hirshfeld_atomic_dipoles = np.array(
        [
            [0.04249319, 0.05486053, 0.00000000],
            [0.13710134, -0.00105126, 0.00000000],
            [-0.03534982, 0.13248706, 0.00000000],
        ]
    )
    assert np.allclose(molecular_calc_chunk.hirshfeld_atomic_dipoles, hirshfeld_atomic_dipoles)
    assert np.allclose(
        molecular_calc_chunk.results["hirshfeld_atomic_dipoles"],
        hirshfeld_atomic_dipoles,
    )


def test_molecular_calc_hirshfeld_dipole(molecular_calc_chunk, molecular_hirshfeld_charges, molecular_positions):
    hirshfeld_dipole = np.sum(molecular_hirshfeld_charges.reshape((-1, 1)) * molecular_positions, axis=1)
    assert np.allclose(molecular_calc_chunk.hirshfeld_dipole, hirshfeld_dipole)
    assert np.allclose(molecular_calc_chunk.results["hirshfeld_dipole"], hirshfeld_dipole)
