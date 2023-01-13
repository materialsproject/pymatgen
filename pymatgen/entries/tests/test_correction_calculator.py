from __future__ import annotations

import os
import unittest

import pytest

from pymatgen.entries.correction_calculator import CorrectionCalculator
from pymatgen.util.testing import PymatgenTest


class CorrectionCalculatorTest(unittest.TestCase):
    def setUp(self):

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

        self.test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "correction_calculator")

    def tearDown(self):
        pass

    def test_normal_corrections(self):
        """
        Test that the values in MPCompatiblity.yaml are reproduced correctly.
        """
        exp_path = os.path.join(self.test_dir, "exp_compounds_norm.json.gz")
        calc_path = os.path.join(self.test_dir, "calc_compounds_norm.json.gz")

        calculator = CorrectionCalculator(exclude_polyanions=self.exclude_polyanions)
        corrs = calculator.compute_from_files(exp_path, calc_path)

        assert corrs == self.normal_corrections

    def test_warnings_options(self):
        """
        Test that compounds can be included/excluded using the allow_{warning} optional parameters.
        """
        exp_path = os.path.join(self.test_dir, "exp_compounds_norm.json.gz")
        calc_path = os.path.join(self.test_dir, "calc_compounds_norm.json.gz")

        calculator = CorrectionCalculator(max_error=1, exclude_polyanions=[], allow_unstable=True)
        corrs = calculator.compute_from_files(exp_path, calc_path)

        assert corrs == self.warnings_allowed_corrections

    def test_no_uncertainties(self):
        """
        Test that corrections can be calculated with no uncertainties.
        """
        exp_path = os.path.join(self.test_dir, "exp_no_error_compounds.json.gz")
        calc_path = os.path.join(self.test_dir, "calc_compounds_norm.json.gz")

        calculator = CorrectionCalculator(exclude_polyanions=self.exclude_polyanions)
        corrs = calculator.compute_from_files(exp_path, calc_path)

        assert corrs == self.no_uncertainties_corrections

    def test_missing_entry_response(self):
        """
        Test that correct error is raised (ValueError) if the input is missing a computed entry.
        """
        exp_path = os.path.join(self.test_dir, "exp_compounds_norm.json.gz")
        calc_path = os.path.join(self.test_dir, "calc_missing_compounds.json.gz")

        calculator = CorrectionCalculator(exclude_polyanions=self.exclude_polyanions)
        with pytest.raises(ValueError):
            calculator.compute_from_files(exp_path, calc_path)


if __name__ == "__main__":
    unittest.main()
