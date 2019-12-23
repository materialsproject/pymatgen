import os
import unittest

from pymatgen.entries.corrections_calc import CorrectionCalculator


class CorrectionCalculatorTest(unittest.TestCase):
    def setUp(self):
        self.normal_corrections = {
            "oxide": (-0.719, 0.0016),
            "peroxide": (-0.47, 0.0187),
            "superoxide": (-0.156, 0.0075),
            "F": (-0.457, 0.0025),
            "Cl": (-0.583, 0.0016),
            "Br": (-0.538, 0.0022),
            "I": (-0.384, 0.0059),
            "N": (-0.332, 0.0084),
            "S": (-0.631, 0.0122),
            "Se": (-0.47, 0.0341),
            "Si": (0.237, 0.0243),
            "Sb": (-0.186, 0.0089),
            "Te": (-0.42, 0.0241),
            "V": (-1.612, 0.0059),
            "Cr": (-1.993, 0.0093),
            "Mn": (-1.61, 0.005),
            "Fe": (-2.252, 0.008),
            "Co": (-1.664, 0.0059),
            "Ni": (-2.732, 0.0105),
            "W": (-4.571, 0.0181),
            "Mo": (-3.075, 0.0084),
            "H": (0.817, 0.0013),
        }

        self.warnings_allowed_corrections = {
            "oxide": (-0.622, 0.0015),
            "peroxide": (-0.471, 0.019),
            "superoxide": (-0.156, 0.0075),
            "F": (-0.454, 0.0025),
            "Cl": (-0.583, 0.0016),
            "Br": (-0.538, 0.0022),
            "I": (-0.383, 0.0059),
            "N": (-0.029, 0.0077),
            "S": (0.529, 0.0075),
            "Se": (-0.319, 0.0304),
            "Si": (0.094, 0.0189),
            "Sb": (-0.19, 0.0089),
            "Te": (-0.421, 0.0238),
            "V": (-1.815, 0.0057),
            "Cr": (-2.204, 0.009),
            "Mn": (-1.732, 0.0049),
            "Fe": (-2.302, 0.0074),
            "Co": (-1.656, 0.0055),
            "Ni": (-2.611, 0.0085),
            "W": (-5.009, 0.0176),
            "Mo": (-3.352, 0.0082),
            "H": (0.808, 0.0013),
        }

        self.no_uncertainties_corrections = {
            "oxide": (-0.732, 0.0338),
            "peroxide": (-0.464, 0.0847),
            "superoxide": (-0.153, 0.0939),
            "F": (-0.403, 0.0349),
            "Cl": (-0.51, 0.0457),
            "Br": (-0.432, 0.0506),
            "I": (-0.31, 0.06),
            "N": (-0.35, 0.1062),
            "S": (-0.638, 0.0664),
            "Se": (-0.358, 0.0939),
            "Si": (-0.366, 0.1716),
            "Sb": (-0.278, 0.1334),
            "Te": (-0.446, 0.1117),
            "V": (-1.668, 0.1312),
            "Cr": (-1.931, 0.1201),
            "Mn": (-1.578, 0.1208),
            "Fe": (-2.203, 0.1071),
            "Co": (-1.697, 0.1366),
            "Ni": (-2.759, 0.1703),
            "W": (-4.343, 0.2046),
            "Mo": (-3.058, 0.1919),
            "H": (1.097, 0.0572),
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
