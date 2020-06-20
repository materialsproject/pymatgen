import os
import unittest

from pymatgen.entries.corrections_calc import CorrectionCalculator


class CorrectionCalculatorTest(unittest.TestCase):
    def setUp(self):
        self.normal_corrections = {
            "oxide": (-0.738, 0.0017),
            "peroxide": (-0.466, 0.0186),
            "superoxide": (-0.16, 0.0075),
            "sulfide": (-0.637, 0.0121),
            "F": (-0.481, 0.0025),
            "Cl": (-0.586, 0.0016),
            "Br": (-0.536, 0.0022),
            "I": (-0.377, 0.0059),
            "N": (-0.344, 0.0084),
            "Se": (-0.473, 0.0341),
            "Si": (-0.406, 0.0621),
            "Sb": (-0.192, 0.0089),
            "Te": (-0.419, 0.024),
            "V": (-1.604, 0.006),
            "Cr": (-1.899, 0.0093),
            "Mn": (-1.592, 0.0051),
            "Fe": (-2.181, 0.009),
            "Co": (-1.521, 0.0059),
            "Ni": (-2.499, 0.0105),
            "W": (-4.584, 0.0181),
            "Mo": (-3.066, 0.0085),
            "H": (-0.178, 0.0013),
        }

        self.warnings_allowed_corrections = {
            "oxide": (-0.633, 0.0015),
            "peroxide": (-0.467, 0.019),
            "superoxide": (-0.16, 0.0075),
            "sulfide": (0.492, 0.0075),
            "F": (-0.478, 0.0025),
            "Cl": (-0.587, 0.0016),
            "Br": (-0.536, 0.0022),
            "I": (-0.376, 0.0059),
            "N": (-0.024, 0.0077),
            "Se": (-0.343, 0.0304),
            "Si": (-0.004, 0.0264),
            "Sb": (-0.197, 0.0089),
            "Te": (-0.42, 0.0238),
            "V": (-1.819, 0.0058),
            "Cr": (-2.128, 0.009),
            "Mn": (-1.722, 0.0049),
            "Fe": (-2.173, 0.0082),
            "Co": (-1.513, 0.0055),
            "Ni": (-2.373, 0.0086),
            "W": (-5.058, 0.0176),
            "Mo": (-3.364, 0.0082),
            "H": (-0.176, 0.0013),
        }

        self.no_uncertainties_corrections = {
            "oxide": (-0.744, 0.0193),
            "peroxide": (-0.466, 0.0478),
            "superoxide": (-0.155, 0.0529),
            "sulfide": (-0.638, 0.0375),
            "F": (-0.436, 0.02),
            "Cl": (-0.517, 0.0258),
            "Br": (-0.436, 0.0285),
            "I": (-0.315, 0.0338),
            "N": (-0.357, 0.0599),
            "Se": (-0.361, 0.053),
            "Si": (-0.439, 0.1017),
            "Sb": (-0.286, 0.0753),
            "Te": (-0.438, 0.063),
            "V": (-1.66, 0.0759),
            "Cr": (-1.85, 0.0679),
            "Mn": (-1.589, 0.0689),
            "Fe": (-2.165, 0.0606),
            "Co": (-1.563, 0.0771),
            "Ni": (-2.516, 0.0961),
            "W": (-4.574, 0.1199),
            "Mo": (-3.05, 0.1085),
            "H": (-0.151, 0.0323),
        }

        self.test_dir = os.path.join(
            os.path.dirname(__file__),
            "..",
            "..",
            "..",
            "test_files",
            "correction_calculator",
        )

    def tearDown(self):
        pass

    def test_normal_corrections(self):
        """
        Test that the values in MPCompatiblity.yaml are reproduced correctly.
        """

        exp_path = os.path.join(self.test_dir, "exp_compounds_norm.gz")
        calc_path = os.path.join(self.test_dir, "calc_compounds_norm.gz")

        calculator = CorrectionCalculator(exp_path, calc_path)
        corrs = calculator.compute_corrections()

        self.assertDictEqual(corrs, self.normal_corrections)

    def test_warnings_options(self):
        """
        Test that compounds can be included/excluded using the allow_{warning} optional parameters.
        """

        exp_path = os.path.join(self.test_dir, "exp_compounds_norm.gz")
        calc_path = os.path.join(self.test_dir, "calc_compounds_norm.gz")

        calculator = CorrectionCalculator(exp_path, calc_path)
        corrs = calculator.compute_corrections(
            allow_polyanions=True, allow_large_errors=True, allow_unstable=True
        )

        self.assertDictEqual(corrs, self.warnings_allowed_corrections)

    def test_no_uncertainties(self):
        """
        Test that corrections can be calculated with no uncertainties.
        """

        exp_path = os.path.join(self.test_dir, "exp_no_error_compounds.gz")
        calc_path = os.path.join(self.test_dir, "calc_compounds_norm.gz")

        calculator = CorrectionCalculator(exp_path, calc_path)
        corrs = calculator.compute_corrections()

        self.assertDictEqual(corrs, self.no_uncertainties_corrections)

    def test_missing_entry_response(self):
        """
        Test that correct error is raised (ValueError) if the input is missing a computed entry.
        """

        exp_path = os.path.join(self.test_dir, "exp_compounds_norm.gz")
        calc_path = os.path.join(self.test_dir, "calc_missing_compounds.gz")

        calculator = CorrectionCalculator(exp_path, calc_path)
        self.assertRaises(ValueError, calculator.compute_corrections)


if __name__ == "__main__":
    unittest.main()
