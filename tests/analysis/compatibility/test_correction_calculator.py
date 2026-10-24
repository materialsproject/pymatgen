from __future__ import annotations

import warnings

import numpy as np
import pytest
from monty.serialization import loadfn

from pymatgen.analysis.compatibility.correction_calculator import CorrectionCalculator
from pymatgen.core import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from tests.testing import TEST_FILES_DIR

TEST_DIR = f"{TEST_FILES_DIR}/entries/correction_calculator"


class TestCorrectionCalculator:
    def setup_method(self):
        self.exclude_polyanions = [
            "SO4",
            "CO3",
            "NO3",
            "OCl3",
            "SiO4",
            "SeO3",
            "TiO3",
            "TiO4",
        ]

        self.normal_corrections = {
            "oxide": (-0.74, 0.0017),
            "peroxide": (-0.466, 0.0186),
            "superoxide": (-0.16, 0.0075),
            "S": (-0.639, 0.0121),
            "F": (-0.485, 0.0025),
            "Cl": (-0.593, 0.0016),
            "Br": (-0.538, 0.0022),
            "I": (-0.377, 0.0059),
            "N": (-0.344, 0.0084),
            "Se": (-0.473, 0.0341),
            "Si": (-0.406, 0.0621),
            "Sb": (-0.192, 0.0089),
            "Te": (-0.419, 0.024),
            "V": (-1.602, 0.006),
            "Cr": (-1.893, 0.0093),
            "Mn": (-1.588, 0.0051),
            "Fe": (-2.182, 0.009),
            "Co": (-1.535, 0.0059),
            "Ni": (-2.504, 0.0105),
            "W": (-4.57, 0.0181),
            "Mo": (-3.058, 0.0085),
            "H": (-0.178, 0.0013),
            "ozonide": (0, 0),
        }

        self.warnings_allowed_corrections = {
            "oxide": (-0.589, 0.0013),
            "peroxide": (-0.467, 0.019),
            "superoxide": (-0.16, 0.0075),
            "S": (-0.637, 0.0123),
            "F": (-0.481, 0.0025),
            "Cl": (-0.594, 0.0016),
            "Br": (-0.538, 0.0022),
            "I": (-0.376, 0.0059),
            "N": (-0.042, 0.0077),
            "Se": (-0.378, 0.0304),
            "Si": (-0.09, 0.0264),
            "Sb": (-0.197, 0.0089),
            "Te": (-0.42, 0.0238),
            "V": (-1.914, 0.0056),
            "Cr": (-2.227, 0.0089),
            "Mn": (-1.781, 0.0048),
            "Fe": (-2.214, 0.0081),
            "Co": (-1.549, 0.0055),
            "Ni": (-2.294, 0.0085),
            "W": (-5.263, 0.0173),
            "Mo": (-3.49, 0.008),
            "H": (-0.176, 0.0013),
            "ozonide": (0, 0),
        }

        self.no_uncertainties_corrections = {
            "oxide": (-0.749, 0.0198),
            "peroxide": (-0.466, 0.0492),
            "superoxide": (-0.155, 0.0545),
            "S": (-0.639, 0.0386),
            "F": (-0.444, 0.0206),
            "Cl": (-0.518, 0.0266),
            "Br": (-0.439, 0.0294),
            "I": (-0.29, 0.0327),
            "N": (-0.357, 0.0617),
            "Se": (-0.361, 0.0546),
            "Si": (-0.439, 0.1048),
            "Sb": (-0.286, 0.0775),
            "Te": (-0.457, 0.0524),
            "V": (-1.655, 0.0782),
            "Cr": (-1.836, 0.0699),
            "Mn": (-1.58, 0.071),
            "Fe": (-2.168, 0.0624),
            "Co": (-1.605, 0.0794),
            "Ni": (-2.516, 0.099),
            "W": (-4.553, 0.1235),
            "Mo": (-3.032, 0.1118),
            "H": (-0.137, 0.0313),
            "ozonide": (0, 0),
        }

    def test_normal_corrections(self):
        """Test that the values in MPCompatibility.yaml are reproduced correctly."""
        exp_path = f"{TEST_DIR}/exp_compounds_norm.json.gz"
        calc_path = f"{TEST_DIR}/calc_compounds_norm.json.gz"

        calculator = CorrectionCalculator(exclude_polyanions=self.exclude_polyanions)
        corrections = calculator.compute_from_files(exp_path, calc_path)

        assert corrections == self.normal_corrections

    def test_warnings_options(self):
        """Test that compounds can be included/excluded using the allow_{warning} optional parameters."""
        exp_path = f"{TEST_DIR}/exp_compounds_norm.json.gz"
        calc_path = f"{TEST_DIR}/calc_compounds_norm.json.gz"

        calculator = CorrectionCalculator(max_error=1, exclude_polyanions=[], allow_unstable=True)
        corrections = calculator.compute_from_files(exp_path, calc_path)

        assert corrections == self.warnings_allowed_corrections

    def test_no_uncertainties(self):
        """Test that corrections can be calculated with no uncertainties."""
        exp_path = f"{TEST_DIR}/exp_no_error_compounds.json.gz"
        calc_path = f"{TEST_DIR}/calc_compounds_norm.json.gz"

        calculator = CorrectionCalculator(exclude_polyanions=self.exclude_polyanions)
        corrections = calculator.compute_from_files(exp_path, calc_path)

        assert corrections == self.no_uncertainties_corrections

    def test_missing_entry_response(self):
        """Test that correct error is raised (ValueError) if the input is missing a computed entry."""
        exp_path = f"{TEST_DIR}/exp_compounds_norm.json.gz"
        calc_path = f"{TEST_DIR}/calc_missing_compounds.json.gz"

        calculator = CorrectionCalculator(exclude_polyanions=self.exclude_polyanions)
        with pytest.raises(ValueError, match="Computed entries missing O"):
            calculator.compute_from_files(exp_path, calc_path)


# Exact corrections (eV per S / F atom) baked into the synthetic dataset below.
SYNTH_CORRECTIONS = {"S": -0.6, "F": -0.4}


def _synthetic_data() -> tuple[list[dict], dict[str, ComputedEntry]]:
    """Build a small dataset whose exp - calc errors are exactly SYNTH_CORRECTIONS.

    Elemental references have zero energy, so the reaction energy equals the compound energy.
    LiS is flagged as a polysulfide and must receive no S correction.
    """
    calc = {
        "Li": ComputedEntry("Li", 0),
        "S": ComputedEntry("S", 0),
        "F2": ComputedEntry("F2", 0),
        "Li2S": ComputedEntry("Li2S", -4.0, data={"e_above_hull": 0}),
        "LiF": ComputedEntry("LiF", -6.0, data={"e_above_hull": 0}),
        "LiS": ComputedEntry("LiS", -2.0, data={"e_above_hull": 0, "sulfide_type": "polysulfide"}),
        "Li3SF": ComputedEntry("Li3SF", -10.0, data={"e_above_hull": 0.05, "sulfide_type": "sulfide"}),
    }
    exp = [
        {"formula": "Li2S", "exp energy": -4.6, "uncertainty": 0.01},
        {"formula": "LiF", "exp energy": -6.4, "uncertainty": 0.02},
        {"formula": "LiS", "exp energy": -2.0, "uncertainty": 0.01},
        {"formula": "Li3SF", "exp energy": -11.0, "uncertainty": 0.03},
    ]
    return exp, calc


@pytest.fixture(scope="module")
def fitted_calculator() -> CorrectionCalculator:
    calculator = CorrectionCalculator(exclude_polyanions="SO4 CO3 NO3 OCl3 SiO4 SeO3 TiO3 TiO4".split())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        calculator.compute_from_files(
            f"{TEST_DIR}/exp_compounds_norm.json.gz", f"{TEST_DIR}/calc_compounds_norm.json.gz"
        )
    return calculator


@pytest.mark.parametrize(("allow_unstable", "expected"), [(False, 0.1), (True, True), (0.5, 0.5)])
def test_allow_unstable_init(allow_unstable, expected):
    assert CorrectionCalculator(allow_unstable=allow_unstable).allow_unstable == expected


def test_synthetic_corrections_exact():
    exp, calc = _synthetic_data()
    calculator = CorrectionCalculator(species=["S", "F"], exclude_polyanions=[])
    corrections = calculator.compute_corrections(exp, calc)

    assert calculator.names == ["Li2S", "LiF", "LiS", "Li3SF"]
    assert calculator.corrections == pytest.approx([SYNTH_CORRECTIONS["S"], SYNTH_CORRECTIONS["F"]], abs=1e-8)
    assert corrections["S"][0] == SYNTH_CORRECTIONS["S"]
    assert corrections["F"][0] == SYNTH_CORRECTIONS["F"]
    assert corrections["ozonide"] == (0, 0)
    # polysulfide LiS gets no S coefficient; Li2S (default) and Li3SF (explicit) are sulfides
    assert calculator.sulfides == ["Li2S", "Li3SF"]
    assert calculator.coeff_mat[2] == [0, 0]
    # residuals vanish for an exactly consistent dataset
    fig = calculator.graph_residual_error()
    assert np.allclose(fig.data[0].y, 0, atol=1e-8)


@pytest.mark.parametrize(
    ("case", "excluded", "match"),
    [
        ("missing", "LiF", "Compound LiF is not found in provided computed entries"),
        ("uncertainty", "LiF", "Compound LiF is excluded from the fit due to high experimental uncertainty"),
        ("polyanion", "Li3SF", "Compound Li3SF contains the poly"),
        ("unstable", "Li3SF", r"Compound Li3SF is unstable and excluded from the fit \(e_above_hull = 0.05\)"),
    ],
)
def test_synthetic_exclusions(case, excluded, match):
    exp, calc = _synthetic_data()
    kwargs: dict = {"species": ["S", "F"], "exclude_polyanions": []}
    if case == "missing":
        del calc[excluded]
    elif case == "uncertainty":
        exp[1]["uncertainty"] = 1.0
    elif case == "polyanion":
        kwargs["exclude_polyanions"] = ["SF"]
    else:
        kwargs["allow_unstable"] = 0.01

    calculator = CorrectionCalculator(**kwargs)
    with pytest.warns(UserWarning, match=match):
        calculator.compute_corrections(exp, calc)

    assert excluded not in calculator.names
    assert len(calculator.names) == 3
    # remaining compounds still determine the corrections exactly
    assert calculator.corrections == pytest.approx([SYNTH_CORRECTIONS["S"], SYNTH_CORRECTIONS["F"]], abs=1e-8)


def test_missing_e_above_hull():
    exp, calc = _synthetic_data()
    calc["LiF"].data.pop("e_above_hull")
    calculator = CorrectionCalculator(species=["S", "F"], exclude_polyanions=[])
    with pytest.raises(ValueError, match="Missing e above hull data"):
        calculator.compute_corrections(exp, calc)


def test_undetectable_specie():
    exp, calc = _synthetic_data()
    calculator = CorrectionCalculator(species=["S", "foo"], exclude_polyanions=[])
    with pytest.raises(ValueError, match="We can't detect this specie='foo' in name='Li2S'"):
        calculator.compute_corrections(exp, calc)


@pytest.mark.parametrize(
    "call",
    [
        lambda calc: calc.graph_residual_error(),
        lambda calc: calc.graph_residual_error_per_species("F"),
        lambda calc: calc.make_yaml(),
    ],
    ids=["graph_residual_error", "graph_residual_error_per_species", "make_yaml"],
)
def test_requires_corrections(call):
    with pytest.raises(RuntimeError, match="Please call compute_corrections or compute_from_files"):
        call(CorrectionCalculator())


def test_graph_residual_error(fitted_calculator):
    calc = fitted_calculator
    fig = calc.graph_residual_error()
    expected = np.abs(np.array(calc.diffs) - np.dot(calc.coeff_mat, calc.corrections))

    assert fig.layout.title.text == "Residual Errors"
    assert fig.layout.yaxis.title.text == "Residual Error (eV/atom)"
    trace = fig.data[0]
    assert len(trace.y) == len(calc.names) == 205
    assert list(trace.x) == list(range(1, 206))
    assert np.allclose(trace.y, np.sort(expected))
    assert sorted(trace.text) == sorted(calc.names)


@pytest.mark.parametrize("specie", ["oxide", "peroxide", "superoxide", "S", "Fe", "H"])
def test_graph_residual_error_per_species(fitted_calculator, specie):
    calc = fitted_calculator
    fig = calc.graph_residual_error_per_species(specie)

    groups = {"oxide": calc.oxides, "peroxide": calc.peroxides, "superoxide": calc.superoxides, "S": calc.sulfides}
    expected = groups.get(specie) or [name for name in calc.names if Composition(name)[specie]]
    assert fig.layout.title.text == f"Residual Errors for {specie}"
    trace = fig.data[0]
    assert sorted(trace.text) == sorted(expected)
    assert len(trace.y) == len(expected) > 0
    assert list(trace.y) == sorted(trace.y)


def test_graph_residual_error_per_species_superoxide(fitted_calculator):
    """Regression: 'superoxide' used to fall through to the sulfide list."""
    fig = fitted_calculator.graph_residual_error_per_species("superoxide")
    assert sorted(fig.data[0].text) == ["CsO2", "KO2", "NaO2"]


def test_graph_residual_error_per_species_invalid(fitted_calculator):
    with pytest.raises(ValueError, match="not a valid specie"):
        fitted_calculator.graph_residual_error_per_species("Xe")


@pytest.mark.parametrize("use_dir", [True, False])
def test_make_yaml(fitted_calculator, tmp_path, monkeypatch, use_dir):
    calc = fitted_calculator
    if use_dir:
        calc.make_yaml(name="Test", dir=str(tmp_path))
    else:
        monkeypatch.chdir(tmp_path)
        calc.make_yaml(name="Test")

    path = tmp_path / "TestCompatibility.yaml"
    assert "# Energy corrections in eV/atom" in path.read_text(encoding="utf-8")
    contents = loadfn(path)

    assert contents["Name"] == "Test"
    corr, uncer = contents["Corrections"], contents["Uncertainties"]
    ggau = ["V", "Cr", "Mn", "Fe", "Co", "Ni", "W", "Mo"]
    for anion in ("O", "F"):
        assert list(corr["GGAUMixingCorrections"][anion]) == ggau
        for el in ggau:
            assert corr["GGAUMixingCorrections"][anion][el] == calc.corrections_dict[el][0]
            assert uncer["GGAUMixingCorrections"][anion][el] == calc.corrections_dict[el][1]
    comp_species = [sp for sp in [*calc.species, "ozonide"] if sp not in ggau]
    assert list(corr["CompositionCorrections"]) == comp_species
    for sp in comp_species:
        assert corr["CompositionCorrections"][sp] == calc.corrections_dict[sp][0]
        assert uncer["CompositionCorrections"][sp] == calc.corrections_dict[sp][1]
    assert corr["CompositionCorrections"]["oxide"] == -0.74
    assert corr["CompositionCorrections"]["ozonide"] == 0
