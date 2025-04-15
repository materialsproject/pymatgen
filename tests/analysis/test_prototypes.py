from __future__ import annotations

import inspect
import re
from itertools import permutations, product
from shutil import which

import pytest

from pymatgen.analysis.prototypes import (
    WYCKOFF_POSITION_RELAB_DICT,
    AflowPrototypeMatcher,
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
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

try:
    import pyxtal
except ImportError:
    pyxtal = None

TEST_DIR = f"{TEST_FILES_DIR}/analysis/prototypes"


class TestAflowPrototypeMatcher(MatSciTest):
    def test_prototype_matching(self):
        af = AflowPrototypeMatcher()

        struct = self.get_structure("Sn")
        prototype = af.get_prototypes(struct)[0]

        assert prototype["tags"] == {
            "aflow": "A_cF8_227_a",
            "mineral": "diamond",
            "pearson": "cF8",
            "strukturbericht": "A4",
        }

        struct = self.get_structure("CsCl")
        prototype = af.get_prototypes(struct)[0]

        assert prototype["tags"] == {
            "aflow": "AB_cP2_221_b_a",
            "mineral": "",
            "pearson": "cP2",
            "strukturbericht": "B2",
        }

        struct = self.get_structure("Li2O")
        prototype = af.get_prototypes(struct)[0]

        assert prototype["tags"] == {
            "aflow": "AB2_cF12_225_a_c",
            "mineral": "Fluorite",
            "pearson": "cF12",
            "strukturbericht": "C1",
        }


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


@pytest.mark.parametrize(("structure", "expected"), zip(TEST_STRUCTS, TEST_PROTOSTRUCTURES, strict=False))
def test_get_protostructure_label_from_spglib(structure, expected):
    """Check that spglib gives correct protostructure label simple cases."""
    assert get_protostructure_label_from_spglib(structure) == expected


def test_get_protostructure_label_from_spglib_edge_case():
    """Check edge case where the symmetry precision is too low."""
    struct = Structure.from_file(f"{TEST_DIR}/U2Pa4Tc6.json")

    defaults = inspect.signature(get_protostructure_label_from_spglib).parameters

    assert defaults["init_symprec"].default == 0.1

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
@pytest.mark.parametrize(("structure", "expected"), zip(TEST_STRUCTS, TEST_PROTOSTRUCTURES, strict=False))
def test_get_protostructure_label_from_aflow(structure, expected):
    """Check that AFLOW CLI gives correct protostructure label simple cases."""
    assert get_protostructure_label_from_aflow(structure, aflow_executable=which("aflow")) == expected


@pytest.mark.parametrize(("structure", "expected"), zip(TEST_STRUCTS, TEST_PROTOSTRUCTURES, strict=False))
def test_get_protostructure_label_from_moyopy(structure, expected):
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
