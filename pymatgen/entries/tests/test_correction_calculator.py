import os
import unittest

from pymatgen.entries.corrections_calc import CorrectionCalculator


class CorrectionCalculatorTest(unittest.TestCase):
    def setUp(self):
        self.normal_corrections = {
            "oxide": (-0.731, 0.0017),
            "peroxide": (-0.476, 0.0186),
            "superoxide": (-0.162, 0.0075),
            "F": (-0.457, 0.0025),
            "Cl": (-0.583, 0.0016),
            "Br": (-0.538, 0.0022),
            "I": (-0.384, 0.0059),
            "N": (-0.34, 0.0084),
            "S": (-0.633, 0.0121),
            "Se": (-0.473, 0.0341),
            "Si": (0.243, 0.0243),
            "Sb": (-0.188, 0.0089),
            "Te": (-0.418, 0.0239),
            "V": (-1.635, 0.006),
            "Cr": (-1.991, 0.0093),
            "Mn": (-1.603, 0.005),
            "Fe": (-2.245, 0.008),
            "Co": (-1.656, 0.0059),
            "Ni": (-2.734, 0.011),
            "Cu": (0.079, 0.0089),
            "Mo": (-3.057, 0.0085),
            "H": (0.819, 0.0013),
        }

        self.warnings_allowed_corrections = {
            "oxide": (-0.717, 0.0015),
            "peroxide": (-0.477, 0.0189),
            "superoxide": (-0.162, 0.0075),
            "F": (-0.446, 0.0025),
            "Cl": (-0.582, 0.0016),
            "Br": (-0.538, 0.0022),
            "I": (-0.383, 0.0059),
            "N": (0.004, 0.0077),
            "S": (0.855, 0.0074),
            "Se": (-0.192, 0.0303),
            "Si": (0.302, 0.0189),
            "Sb": (-0.192, 0.0089),
            "Te": (-0.42, 0.0237),
            "V": (-1.661, 0.0058),
            "Cr": (-1.976, 0.009),
            "Mn": (-1.62, 0.0049),
            "Fe": (-2.222, 0.0074),
            "Co": (-1.637, 0.0055),
            "Ni": (-3.143, 0.0084),
            "Cu": (0.142, 0.0078),
            "Mo": (-3.095, 0.0082),
            "H": (0.81, 0.0013),
        }
        self.no_uncertainties_corrections = {
            "oxide": (-0.744, 0.0357),
            "peroxide": (-0.471, 0.0859),
            "superoxide": (-0.159, 0.0952),
            "F": (-0.419, 0.0363),
            "Cl": (-0.511, 0.0464),
            "Br": (-0.433, 0.0513),
            "I": (-0.31, 0.0608),
            "N": (-0.355, 0.1077),
            "S": (-0.639, 0.0673),
            "Se": (-0.361, 0.0953),
            "Si": (-0.367, 0.174),
            "Sb": (-0.281, 0.1353),
            "Te": (-0.445, 0.1132),
            "V": (-1.726, 0.1362),
            "Cr": (-1.922, 0.1227),
            "Mn": (-1.564, 0.125),
            "Fe": (-2.196, 0.1101),
            "Co": (-1.666, 0.1416),
            "Ni": (-2.743, 0.1783),
            "Cu": (0.065, 0.1177),
            "Mo": (-3.035, 0.2004),
            "H": (1.098, 0.058),
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
