from __future__ import annotations

import inspect
import json
import math
import re
import subprocess
from itertools import permutations, product
from shutil import which
from types import SimpleNamespace

import pandas as pd
import pytest

import pymatgen.analysis.prototypes as prototypes_module
from pymatgen.analysis.prototypes import (
    WYCKOFF_POSITION_RELAB_DICT,
    AflowPrototypeMatcher,
    PrototypeDatabaseMatcher,
    _find_translations,
    count_crystal_dof,
    count_crystal_sites,
    count_distinct_wyckoff_letters,
    count_wyckoff_positions,
    get_anonymous_formula_from_prototype_formula,
    get_formula_from_protostructure_label,
    get_protostructure_label,
    get_protostructure_label_from_aflow,
    get_protostructure_label_from_moyopy,
    get_protostructure_label_from_spg_analyzer,
    get_protostructure_label_from_spglib,
    get_protostructures_from_aflow_label_and_composition,
    get_prototype_formula_from_composition,
    get_prototype_from_protostructure,
    get_random_structure_for_protostructure,
)
from pymatgen.core.structure import Composition, Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

try:
    import pyxtal
except ImportError:
    pyxtal = None

TEST_DIR = f"{TEST_FILES_DIR}/analysis/prototypes"


def assert_common_prototype_matches(test, af):
    struct = test.get_structure("Sn")
    prototype = af.get_prototypes(struct)[0]

    assert prototype["tags"] == {
        "aflow": "A_cF8_227_a",
        "mineral": "diamond",
        "pearson": "cF8",
        "strukturbericht": "A4",
    }

    struct = test.get_structure("CsCl")
    prototype = af.get_prototypes(struct)[0]

    assert prototype["tags"] == {
        "aflow": "AB_cP2_221_b_a",
        "mineral": "",
        "pearson": "cP2",
        "strukturbericht": "B2",
    }

    struct = test.get_structure("Li2O")
    prototype = af.get_prototypes(struct)[0]

    assert prototype["tags"] == {
        "aflow": "AB2_cF12_225_a_c",
        "mineral": "Fluorite",
        "pearson": "cF12",
        "strukturbericht": "C1",
    }


class TestAflowPrototypeMatcher(MatSciTest):
    def test_prototype_matching(self):
        with pytest.warns(DeprecationWarning, match="AflowPrototypeMatcher is deprecated"):
            af = AflowPrototypeMatcher()

        assert_common_prototype_matches(self, af)

    def test_deprecated_aflow_prototype_matcher(self):
        with pytest.warns(DeprecationWarning, match="AflowPrototypeMatcher is deprecated"):
            AflowPrototypeMatcher()

    def test_deprecated_aflow_prototype_matcher_rejects_prototype_db(self):
        with (
            pytest.warns(DeprecationWarning, match="AflowPrototypeMatcher is deprecated"),
            pytest.raises(TypeError, match="unexpected keyword argument 'prototype_db'"),
        ):
            AflowPrototypeMatcher(prototype_db=[])


class TestPrototypeDatabaseMatcher(MatSciTest):
    def test_prototype_matching(self):
        from pymatgen.analysis.prototypes._data import AFLOW_PROTOTYPE_LIBRARY

        assert_common_prototype_matches(self, PrototypeDatabaseMatcher(pd.DataFrame(AFLOW_PROTOTYPE_LIBRARY)))


class TestCustomPrototypeDatabaseMatcher(MatSciTest):
    class CustomPrototypeMatcher(PrototypeDatabaseMatcher):
        @staticmethod
        def _get_structure(row):
            structure = row["structure"]
            return Structure.from_dict(structure) if isinstance(structure, dict) else structure

        @staticmethod
        def _get_entry_data(row):
            structure = row["structure"]
            return {
                "type": row["mineral"],
                "distance": row["distance"],
                "structure": Structure.from_dict(structure) if isinstance(structure, dict) else structure,
            }

    def test_prototype_db_matching(self):
        struct = self.get_structure("Sn")
        prototype_db = pd.DataFrame([{"mineral": "tin", "distance": 0.1, "structure": struct}])
        prototype = self.CustomPrototypeMatcher(prototype_db).get_prototypes(struct)[0]

        assert prototype == {"type": "tin", "distance": 0.1, "structure": struct}

    def test_prototype_db_matching_with_structure_dicts(self):
        struct = self.get_structure("Sn")
        prototype_db = pd.DataFrame([{"mineral": "tin", "distance": 0.1, "structure": struct.as_dict()}])
        prototype = self.CustomPrototypeMatcher(prototype_db).get_prototypes(struct)[0]

        assert prototype == {"type": "tin", "distance": 0.1, "structure": struct}


PROTOSTRUCTURE_SET = [
    ("A20BC14D8E5F2_oP800_61_40c_2c_28c_16c_10c_4c:C-Cd-H-N-O-S"),
    ("ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si"),
    ("ABC6D2_mC40_15_e_a_3f_f:Ca-Fe-O-Si"),
    ("A6B11CD7_aP50_2_6i_ac10i_i_7i:C-H-N-O"),
    ("ABC2D2_mC24_15_e_e_f_f:Ca-Fe-O-Si"),
    ("A3B2CD4_tP10_115_ag_g_b_cdg:Al-Ce-Ga-Pd"),
    ("AB2_cF576_228_h_fgh:Ba-Ti"),
    ("AB3C_cP5_221_a_c_b:Ba-O-Ti"),
]

TEST_STRUCTS = [
    Structure(  # NaCl structure
        lattice=[[2, 2, 0], [0, 2, 2], [2, 0, 2]],
        species=["Na", "Cl"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    ),
    Structure(  # CsCl structure
        lattice=[[4, 0, 0], [0, 4, 0], [0, 0, 4]],
        species=["Cs", "Cl"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    ),
    Structure(  # ZnO zincblende structure
        lattice=[[2, 2, 0], [0, 2, 2], [2, 0, 2]],
        species=["Zn", "O"],
        coords=[[0, 0, 0], [0.25, 0.25, 0.25]],
    ),
    Structure(  # ZnO wurtzite structure
        lattice=Lattice.from_parameters(a=3.8227, b=3.8227, c=6.2607, alpha=90, beta=90, gamma=120),
        species=["Zn", "O", "Zn", "O"],
        coords=[
            [1 / 3, 2 / 3, 0],
            [2 / 3, 1 / 3, 0.3748],
            [2 / 3, 1 / 3, 1 / 2],
            [1 / 3, 2 / 3, 1 / 2 + 0.3748],
        ],
    ),
    Structure(
        lattice=[[3.9, 0, 0], [0, 3.9, 0], [0, 0, 3.9]],
        species=["Sr", "Ti", "O", "O", "O"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
    ),
    Structure(
        lattice=[[5.76, 0, 0], [0, 5.76, 0], [0, 0, 5.76]],
        species=["Al", "Fe", "Fe", "Fe", "Al", "Fe", "Fe", "Fe"],
        coords=[
            [0, 0, 0],
            [0.25, 0.25, 0.25],
            [0.5, 0.5, 0],
            [0.75, 0.75, 0.25],
            [0, 0.5, 0.5],
            [0.25, 0.75, 0.75],
            [0.5, 0, 0.5],
            [0.75, 0.25, 0.75],
        ],
    ),
    Structure.from_file(f"{TEST_DIR}/ABC6D2_mC40_15_e_e_3f_f.cif"),
]

TEST_PROTOSTRUCTURES = [
    "AB_cF8_225_a_b:Cl-Na",
    "AB_cP2_221_a_b:Cl-Cs",
    "AB_cF8_216_a_c:O-Zn",
    "AB_hP4_186_b_b:O-Zn",
    "A3BC_cP5_221_c_a_b:O-Sr-Ti",
    "AB3_tP4_115_a_cg:Al-Fe",
    "ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si",
]


@pytest.mark.parametrize(("structure", "expected"), list(zip(TEST_STRUCTS, TEST_PROTOSTRUCTURES, strict=True)))
def test_get_protostructure_label_from_spglib(structure, expected):
    """Check that spglib gives correct protostructure label simple cases."""
    assert get_protostructure_label_from_spglib(structure) == expected


def test_get_protostructure_label_from_spglib_edge_case():
    """Check edge case where the symmetry precision is too low."""
    struct = Structure.from_file(f"{TEST_DIR}/U2Pa4Tc6.json")

    defaults = inspect.signature(get_protostructure_label_from_spglib).parameters

    assert math.isclose(defaults["init_symprec"].default, 0.1)

    spg_analyzer = SpacegroupAnalyzer(struct, symprec=defaults["init_symprec"].default, angle_tolerance=5)

    raises_str = "Invalid WP multiplicities - A2B3C_hP6_191_c_2g_a:Pa-Tc-U, expected U(PaTc3)2 to be UPa2Tc3"
    with pytest.raises(ValueError, match=re.escape(raises_str)):
        get_protostructure_label_from_spg_analyzer(spg_analyzer, raise_errors=True)

    assert get_protostructure_label_from_spg_analyzer(spg_analyzer, raise_errors=False) == raises_str

    # Test that it gives invalid protostructure if fallback is None.
    with pytest.raises(ValueError, match=re.escape(raises_str)):
        get_protostructure_label_from_spglib(struct, raise_errors=True, fallback_symprec=None)

    assert get_protostructure_label_from_spglib(struct, raise_errors=False, fallback_symprec=None) == raises_str

    assert get_protostructure_label_from_spglib(struct, raise_errors=True) == ("A2B3C_hP6_191_c_g_a:Pa-Tc-U")

    assert get_protostructure_label_from_spglib(struct, raise_errors=False) == ("A2B3C_hP6_191_c_g_a:Pa-Tc-U")


@pytest.mark.parametrize(
    ("protostructure_label", "expected"),
    [
        ("ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si", 6),  # esseneite
        ("A6B11CD7_aP50_2_6i_ac10i_i_7i:C-H-N-O", 26),
        ("foo_bar_47_abc_A_b:X-Y-Z", 5),
    ],
)
def test_count_wyckoff_positions(protostructure_label, expected):
    count = count_wyckoff_positions(protostructure_label)
    assert isinstance(count, int)
    assert count == expected


def test_count_crystal_dof():
    """Count the number of coarse-grained parameters in esseneite"""
    count = count_crystal_dof("ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si")
    assert isinstance(count, int)
    assert count == 18


def test_count_crystal_sites():
    """Count the number of sites in esseneite"""
    count = count_crystal_sites("ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si")
    assert isinstance(count, int)
    assert count == 40


@pytest.mark.parametrize(
    ("protostructure_label", "expected"),
    [
        ("ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si", "ABC2D6_mC40_15_e_e_f_3f"),
        ("ABC6D2_mC40_15_e_a_3f_f:Ca-Fe-O-Si", "ABC2D6_mC40_15_a_e_f_3f"),
        # this case failing means doesn't handle single element materials
        ("A_tI8_141_ea:Ca", "A_tI8_141_ae"),  # not ea is non-canonical
        # this case failing means not reordering elements based on int but first digit
        ("A4BC20D2_oC108_41_2b_a_10b_b:B-Ca-H-N", "AB2C4D20_oC108_41_a_b_2b_10b"),
    ],
)
def test_get_prototype_from_protostructure(protostructure_label, expected):
    """Get a recanonicalized prototype string without chemical system"""
    aflow_label, chemsys = protostructure_label.split(":")
    prototype_formula, pearson_symbol, spg_num, *wyckoffs = aflow_label.split("_")

    element_wyckoff = "_".join(wyckoffs)

    isopointal_element_wyckoffs = list(
        {element_wyckoff.translate(str.maketrans(trans)) for trans in WYCKOFF_POSITION_RELAB_DICT[spg_num]}
    )

    protostructure_labels = [
        f"{prototype_formula}_{pearson_symbol}_{spg_num}_{element_wyckoff}:{chemsys}"
        for element_wyckoff in isopointal_element_wyckoffs
    ]

    assert all(
        get_prototype_from_protostructure(protostructure_label) == expected
        for protostructure_label in protostructure_labels
    )


@pytest.mark.parametrize(
    ("aflow_label", "composition", "expected"),
    [
        (
            "AB2C3D4_tP10_115_a_g_bg_cdg",
            "Ce2Al3GaPd4",
            "A3B2CD4_tP10_115_ag_g_b_cdg:Al-Ce-Ga-Pd",
        ),
        # checks that we can handle cases where one element could be on multiple sites
        (
            "ABC3_oP20_62_a_c_cd",
            "YbNiO3",
            "AB3C_oP20_62_c_cd_a:Ni-O-Yb AB3C_oP20_62_a_cd_c:Ni-O-Yb",
        ),
    ],
)
def test_get_protostructures_from_aflow_label_and_composition(aflow_label, composition, expected):
    protostructures = get_protostructures_from_aflow_label_and_composition(aflow_label, Composition(composition))
    assert set(protostructures) == set(expected.split(" "))

    # check the round trip
    assert all(get_prototype_from_protostructure(protostructure) == aflow_label for protostructure in protostructures)


@pytest.mark.parametrize(
    ("dict1", "dict2", "expected"),
    [
        # Test case 1: Simple valid translation
        ({"a": 1, "b": 2}, {"x": 1, "y": 2}, [{"a": "x", "b": "y"}]),
        # Test case 2: Multiple valid translations
        (
            {"a": 1, "b": 1, "c": 1},
            {"x": 1, "y": 1, "z": 1},
            [dict(zip(["a", "b", "c"], perm, strict=False)) for perm in permutations(["x", "y", "z"])],
        ),
        # Test case 3: No valid translations (different values)
        ({"a": 1, "b": 2}, {"x": 1, "y": 3}, []),
        # Test case 4: No valid translations (different number of items)
        ({"a": 1, "b": 2}, {"x": 1, "y": 2, "z": 3}, []),
        # Test case 5: Empty dictionaries
        ({}, {}, [{}]),
        # Test case 6: Larger dictionaries
        (
            {"a": 1, "b": 4, "c": 3, "d": 4},
            {"w": 4, "x": 3, "y": 4, "z": 1},
            [
                {"a": "z", "b": "y", "c": "x", "d": "w"},
                {"a": "z", "b": "w", "c": "x", "d": "y"},
            ],
        ),
    ],
)
def test_find_translations(dict1, dict2, expected):
    result = _find_translations(dict1, dict2)
    assert len(result) == len(expected)
    for translation in result:
        assert translation in expected


# Additional test for performance with larger input
def test_find_translations_performance():
    dict1 = {f"key{idx}": idx for idx in range(8)}
    dict2 = {f"val{idx}": idx for idx in range(8)}
    result = _find_translations(dict1, dict2)
    assert len(result) == 1  # There should be only one valid translation


@pytest.mark.parametrize(
    ("composition", "expected"),
    [("Ce2Al3GaPd4", "A3B2CD4"), ("YbNiO3", "AB3C"), ("K2NaAlF6", "AB6C2D")],
)
def test_get_prototype_formula_from_composition(composition: str, expected: str):
    assert get_prototype_formula_from_composition(Composition(composition)) == expected


@pytest.mark.parametrize(
    ("protostructure_label", "expected"),
    [("AB3C_oP20_62_c_cd_a:Ni-O-Yb", "NiO3Yb")],
)
def test_get_formula_from_protostructure_label(protostructure_label: str, expected: str):
    assert get_formula_from_protostructure_label(protostructure_label) == expected


@pytest.mark.parametrize(
    ("anonymous_formula", "prototype_formula"),
    [("AB", "AB"), ("A2B", "AB2"), ("A3B2CD4", "AB2C3D4")],
)
def test_get_anonymous_formula_from_prototype_formula(anonymous_formula: str, prototype_formula: str):
    assert get_anonymous_formula_from_prototype_formula(anonymous_formula) == prototype_formula


@pytest.mark.parametrize(
    ("protostructure_label", "expected"),
    [
        ("A20BC14D8E5F2_oP800_61_40c_2c_28c_16c_10c_4c:C-Cd-H-N-O-S", 1),
        ("ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si", 2),
        ("ABC6D2_mC40_15_e_a_3f_f:Ca-Fe-O-Si", 3),
        ("A6B11CD7_aP50_2_6i_ac10i_i_7i:C-H-N-O", 3),
    ],
)
def test_count_distinct_wyckoff_letters(protostructure_label, expected):
    assert count_distinct_wyckoff_letters(protostructure_label) == expected


@pytest.mark.skipif(which("aflow") is None, reason="AFLOW CLI not installed")
@pytest.mark.parametrize(("structure", "expected"), list(zip(TEST_STRUCTS, TEST_PROTOSTRUCTURES, strict=True)))
def test_get_protostructure_label_from_aflow(structure, expected):
    """Check that AFLOW CLI gives correct protostructure label simple cases."""
    assert get_protostructure_label_from_aflow(structure, aflow_executable=which("aflow")) == expected


@pytest.fixture
def mock_aflow(monkeypatch):
    """Replace the AFLOW CLI with a stub reporting a given AFLOW label.

    Returns a function that installs the stub for ``aflow_label`` and returns the
    list in which each ``subprocess.run`` call is recorded as ``(cmd, kwargs)``.
    """
    calls: list[tuple[list[str], dict]] = []

    def install(aflow_label: str) -> list[tuple[list[str], dict]]:
        def fake_run(cmd, **kwargs):
            calls.append((cmd, kwargs))
            stdout = json.dumps({"aflow_prototype_label": aflow_label})
            return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")

        monkeypatch.setattr(prototypes_module, "which", lambda cmd: cmd or None)
        monkeypatch.setattr(prototypes_module, "subprocess", SimpleNamespace(run=fake_run))
        return calls

    return install


@pytest.mark.parametrize(
    ("structure", "aflow_label", "expected"),
    [
        (TEST_STRUCTS[0], "AB_cF8_225_a_b", "AB_cF8_225_a_b:Cl-Na"),  # rocksalt
        # equivalent origin choice must be canonicalized
        (TEST_STRUCTS[0], "AB_cF8_225_b_a", "AB_cF8_225_a_b:Cl-Na"),
        (TEST_STRUCTS[4], "A3BC_cP5_221_c_a_b", "A3BC_cP5_221_c_a_b:O-Sr-Ti"),  # perovskite
        (TEST_STRUCTS[4], "A3BC_cP5_221_d_b_a", "A3BC_cP5_221_c_a_b:O-Sr-Ti"),
        (TEST_STRUCTS[5], "AB3_tP4_115_a_cg", "AB3_tP4_115_a_cg:Al-Fe"),
        (TEST_STRUCTS[6], "ABC6D2_mC40_15_e_e_3f_f", "ABC6D2_mC40_15_e_e_3f_f:Ca-Fe-O-Si"),  # esseneite
    ],
)
def test_get_protostructure_label_from_aflow_mocked(mock_aflow, structure, aflow_label, expected):
    """Parse (mocked) AFLOW CLI output without needing the aflow binary."""
    calls = mock_aflow(aflow_label)

    assert get_protostructure_label_from_aflow(structure) == expected

    (cmd, kwargs), *_ = calls
    assert len(calls) == 1
    assert cmd == ["aflow", "--prototype", "--print=json", "cat"]
    assert kwargs["input"] == structure.to(fmt="poscar")
    assert kwargs["check"] is True


def test_get_protostructure_label_from_aflow_executable(mock_aflow):
    calls = mock_aflow("AB_cF8_225_a_b")
    label = get_protostructure_label(TEST_STRUCTS[0], method="aflow", aflow_executable="my-aflow")

    assert label == "AB_cF8_225_a_b:Cl-Na"
    assert calls[0][0][0] == "my-aflow"


def test_get_protostructure_label_from_aflow_invalid_multiplicities(mock_aflow):
    # Na on 8c and Cl on 4a of Fm-3m is inconsistent with NaCl
    mock_aflow("AB_cF8_225_a_c")
    err_msg = "Invalid WP multiplicities - AB_cF8_225_a_c:Cl-Na, expected Na2Cl to be NaCl"

    assert get_protostructure_label_from_aflow(TEST_STRUCTS[0]) == err_msg
    with pytest.raises(ValueError, match=re.escape(err_msg)):
        get_protostructure_label_from_aflow(TEST_STRUCTS[0], raise_errors=True)


def test_get_protostructure_label_from_aflow_missing_binary(monkeypatch, tmp_path):
    monkeypatch.setattr(prototypes_module, "which", lambda _cmd: None)
    with pytest.raises(FileNotFoundError, match="AFLOW could not be found"):
        get_protostructure_label_from_aflow(TEST_STRUCTS[0])

    monkeypatch.undo()
    missing = str(tmp_path / "no-such-aflow")
    with pytest.raises(FileNotFoundError, match="AFLOW could not be found"):
        get_protostructure_label_from_aflow(TEST_STRUCTS[0], aflow_executable=missing)


@pytest.mark.parametrize(("structure", "expected"), list(zip(TEST_STRUCTS, TEST_PROTOSTRUCTURES, strict=True)))
def test_get_protostructure_label_from_moyopy(structure, expected):
    pytest.importorskip("moyopy", reason="moyopy not installed")

    """Check that moyopy gives correct protostructure label simple cases."""
    assert get_protostructure_label_from_moyopy(structure) == expected, (
        f"unexpected moyopy protostructure for {structure=}"
    )


@pytest.mark.skipif(pyxtal is None, reason="pyxtal not installed")
@pytest.mark.xfail(reason="pyxtal is non-deterministic and symmetry can increase in random crystal")
@pytest.mark.parametrize(
    "protostructure",
    PROTOSTRUCTURE_SET,
)
def test_moyopy_spglib_consistency(protostructure):
    """Check that moyopy and spglib give consistent results."""
    struct = get_random_structure_for_protostructure(protostructure)

    moyopy_label = get_protostructure_label_from_moyopy(struct)
    spglib_label = get_protostructure_label_from_spglib(struct)

    assert moyopy_label == spglib_label, f"spglib moyopy protostructure mismatch for {protostructure}"


@pytest.mark.skipif(pyxtal is None, reason="pyxtal not installed")
@pytest.mark.xfail(reason="pyxtal is non-deterministic and symmetry can increase in random crystal")
@pytest.mark.parametrize(
    ("protostructure", "method"),
    list(product(PROTOSTRUCTURE_SET, ["spglib", "moyopy"])),
)
def test_get_random_structure_for_protostructure_roundtrip(protostructure: str, method: str):
    """Check roundtrip for generating a random structure from a prototype string"""
    assert protostructure == get_protostructure_label(
        get_random_structure_for_protostructure(protostructure),
        method=method,
    )


@pytest.mark.skipif(pyxtal is None, reason="pyxtal not installed")
@pytest.mark.parametrize(
    "protostructure",
    PROTOSTRUCTURE_SET,
)
def test_get_random_structure_for_protostructure_random(protostructure):
    """Check roundtrip for generating a random structure from a prototype string"""
    s1 = get_random_structure_for_protostructure(protostructure)
    s2 = get_random_structure_for_protostructure(protostructure)

    assert s1.composition == s2.composition
    assert s1.lattice != s2.lattice


SPINEL = Structure.from_spacegroup(  # MgAl2O4 normal spinel
    "Fd-3m",
    Lattice.cubic(8.08),
    ["Al", "Mg", "O"],
    [[0.125, 0.125, 0.125], [0.5, 0.5, 0.5], [0.3875, 0.3875, 0.3875]],
)
RUTILE = Structure.from_spacegroup(  # TiO2
    "P4_2/mnm",
    Lattice.tetragonal(4.594, 2.959),
    ["Ti", "O"],
    [[0, 0, 0], [0.3048, 0.3048, 0]],
)
ROCKSALT, PEROVSKITE = TEST_STRUCTS[0], TEST_STRUCTS[4]

# (structure, protostructure label, prototype, n_sites, n_dof, formula, mineral tag in AFLOW library)
KNOWN_PROTOTYPES = [
    (ROCKSALT, "AB_cF8_225_a_b:Cl-Na", "AB_cF8_225_a_b", 8, 1, "ClNa", "Halite, Rock Salt"),
    (PEROVSKITE, "A3BC_cP5_221_c_a_b:O-Sr-Ti", "ABC3_cP5_221_a_b_c", 5, 1, "O3SrTi", "(Cubic) Perovskite"),
    (SPINEL, "A2BC4_cF56_227_c_b_e:Al-Mg-O", "AB2C4_cF56_227_a_d_e", 56, 2, "Al2MgO4", "Spinel"),
    (RUTILE, "A2B_tP6_136_f_a:O-Ti", "AB2_tP6_136_a_f", 6, 3, "O2Ti", "Rutile"),
]
KNOWN_PROTOTYPE_IDS = ["rocksalt", "perovskite", "spinel", "rutile"]


@pytest.mark.parametrize("method", ["spglib", "moyopy"])
@pytest.mark.parametrize(
    ("structure", "label", "prototype", "n_sites", "n_dof", "formula", "mineral"),
    KNOWN_PROTOTYPES,
    ids=KNOWN_PROTOTYPE_IDS,
)
def test_known_prototype_labels(method, structure, label, prototype, n_sites, n_dof, formula, mineral):
    """Protostructure labels and derived quantities for textbook prototypes."""
    if method == "moyopy":
        pytest.importorskip("moyopy", reason="moyopy not installed")

    assert get_protostructure_label(structure, method=method) == label
    assert get_prototype_from_protostructure(label) == prototype
    assert count_crystal_sites(label) == n_sites
    assert count_crystal_dof(label) == n_dof
    assert get_formula_from_protostructure_label(label) == formula
    assert Composition(formula).reduced_formula == structure.composition.reduced_formula


@pytest.mark.parametrize(
    ("structure", "label", "prototype", "n_sites", "n_dof", "formula", "mineral"),
    KNOWN_PROTOTYPES,
    ids=KNOWN_PROTOTYPE_IDS,
)
def test_protostructures_from_aflow_label_round_trip(structure, label, prototype, n_sites, n_dof, formula, mineral):
    """The anonymous prototype plus the (unreduced) cell composition regenerates the
    protostructure label among the possible decorations.
    """
    protostructures = get_protostructures_from_aflow_label_and_composition(prototype, structure.composition)

    assert label in protostructures
    assert all(get_prototype_from_protostructure(proto) == prototype for proto in protostructures)
    if structure is PEROVSKITE:  # Sr and Ti both 1 per f.u., so either may take the 1a/1b site
        assert sorted(protostructures) == ["A3BC_cP5_221_c_a_b:O-Sr-Ti", "A3BC_cP5_221_c_b_a:O-Sr-Ti"]
    else:
        assert protostructures == [label]


@pytest.fixture(scope="module")
def aflow_library_matcher():
    from pymatgen.analysis.prototypes._data import AFLOW_PROTOTYPE_LIBRARY

    return PrototypeDatabaseMatcher(pd.DataFrame(AFLOW_PROTOTYPE_LIBRARY))


@pytest.mark.parametrize(
    ("structure", "label", "prototype", "n_sites", "n_dof", "formula", "mineral"),
    KNOWN_PROTOTYPES,
    ids=KNOWN_PROTOTYPE_IDS,
)
def test_prototype_database_matcher_known_prototypes(
    aflow_library_matcher, structure, label, prototype, n_sites, n_dof, formula, mineral
):
    matches = aflow_library_matcher.get_prototypes(structure)

    assert [match["tags"]["mineral"] for match in matches] == [mineral]
    # the library's AFLOW label canonicalizes to the same prototype as ours
    aflow_label = matches[0]["tags"]["aflow"]
    assert get_prototype_from_protostructure(f"{aflow_label}:{label.split(':')[1]}") == prototype


def test_prototype_database_matcher_no_match(aflow_library_matcher):
    """A low-symmetry arrangement absent from the library returns None, not []."""
    struct = Structure(
        Lattice.cubic(5),
        ["H", "H", "H", "Li"],
        [[0, 0, 0], [0.1, 0.2, 0.3], [0.5, 0.6, 0.1], [0.7, 0.2, 0.9]],
    )
    assert aflow_library_matcher.get_prototypes(struct) is None


def test_prototype_database_matcher_unresolvable_duplicates():
    """Identical prototypes can't be separated by tightening tolerances; the loop must
    terminate once ltol drops below 0.01 and return all matches.
    """
    prototype_db = pd.DataFrame(
        [
            {"snl": SimpleNamespace(structure=ROCKSALT), "name": "rocksalt-1"},
            {"snl": SimpleNamespace(structure=ROCKSALT), "name": "rocksalt-2"},
        ]
    )
    matches = PrototypeDatabaseMatcher(prototype_db).get_prototypes(ROCKSALT)

    assert [match["name"] for match in matches] == ["rocksalt-1", "rocksalt-2"]


def test_get_protostructure_label_invalid_method():
    with pytest.raises(ValueError, match="Invalid method: foo"):
        get_protostructure_label(ROCKSALT, method="foo")


# spglib's own deprecation of OLD_ERROR_HANDLING, triggered via SpacegroupAnalyzer on failure
@pytest.mark.filterwarnings("ignore:Set OLD_ERROR_HANDLING to false:DeprecationWarning")
def test_get_protostructure_label_from_spglib_undetermined_symmetry(monkeypatch):
    """Coincident atoms make spglib fail; the error is annotated unless raise_errors."""
    # ase's spglib wrapper sets SPGLIB_OLD_ERROR_HANDLING=false and does not unset it, so after any
    # test that goes through ase symmetry code spglib raises its own message instead of returning None.
    monkeypatch.delenv("SPGLIB_OLD_ERROR_HANDLING", raising=False)
    struct = Structure(Lattice.cubic(4), ["Na", "Na"], [[0, 0, 0], [0, 0, 0]], validate_proximity=False)

    assert get_protostructure_label_from_spglib(struct) == "Unable to determine symmetry"
    with pytest.raises(ValueError, match="Unable to determine symmetry"):
        get_protostructure_label_from_spglib(struct, raise_errors=True)


def test_get_protostructure_label_from_moyopy_invalid_multiplicities(monkeypatch):
    pytest.importorskip("moyopy", reason="moyopy not installed")

    # pretend moyopy put Na on 4a and Cl on 8c of Fm-3m, inconsistent with NaCl
    fake_dataset = SimpleNamespace(number=225, pearson_symbol="cF8", orbits=[0, 1], wyckoffs=["4a", "8c"])
    monkeypatch.setattr(prototypes_module, "moyopy", SimpleNamespace(MoyoDataset=lambda _cell, symprec: fake_dataset))
    err_msg = "Invalid WP multiplicities - AB_cF8_225_c_a:Cl-Na, expected NaCl2 to be NaCl"

    assert get_protostructure_label_from_moyopy(ROCKSALT) == err_msg
    with pytest.raises(ValueError, match=re.escape(err_msg)):
        get_protostructure_label_from_moyopy(ROCKSALT, raise_errors=True)


def test_missing_optional_dependencies(monkeypatch):
    monkeypatch.setattr(prototypes_module, "HAS_MOYOPY", False)
    monkeypatch.setattr(prototypes_module, "HAS_PYXTAL", False)

    with pytest.raises(ImportError, match="moyopy not found"):
        get_protostructure_label(ROCKSALT, method="moyopy")
    with pytest.raises(ImportError, match="pyxtal is required"):
        get_random_structure_for_protostructure("AB_cF8_225_a_b:Cl-Na")


@pytest.mark.parametrize(
    ("composition", "expected"),
    [
        ("Na2Cl2", "AB"),  # reduced by gcd
        ("Fe4O6", "A2B3"),
        ("Fe0.5Ni0.5", "A0.5B0.5"),  # non-integer amounts are not reduced
        ("Li0.5CoO2", "AB0.5C2"),
    ],
)
def test_get_prototype_formula_from_composition_reduction(composition, expected):
    assert get_prototype_formula_from_composition(Composition(composition)) == expected
