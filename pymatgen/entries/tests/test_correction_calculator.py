import os
import unittest

from pymatgen.entries.corrections_calc import CorrectionCalculator


class CorrectionCalculatorTest(unittest.TestCase):
    def setUp(self):
        self.normal_corrections = {
            'oxide': (-0.723, 0.0017),
            'peroxide': (-0.469, 0.0186),
            'superoxide': (-0.156, 0.0075),
            'sulfide': (-0.632, 0.0121),
            'F': (-0.46, 0.0025),
            'Cl': (-0.583, 0.0016),
            'Br': (-0.538, 0.0022),
            'I': (-0.384, 0.0059),
            'N': (-0.332, 0.0084),
            'Se': (-0.47, 0.0341),
            'Si': (-0.414, 0.0621),
            'Sb': (-0.186, 0.0089),
            'Te': (-0.421, 0.024),
            'V': (-1.598, 0.006),
            'Cr': (-1.957, 0.0093),
            'Mn': (-1.59, 0.0051),
            'Fe': (-2.232, 0.009),
            'Co': (-1.519, 0.0059),
            'Ni': (-2.583, 0.0105),
            'W': (-4.499, 0.0181),
            'Mo': (-3.063, 0.0085),
            'H': (-0.184, 0.0013)
        }

        self.warnings_allowed_corrections = {
            'oxide': (-0.616, 0.0015),
            'peroxide': (-0.47, 0.019),
            'superoxide': (-0.156, 0.0075),
            'sulfide': (0.534, 0.0075),
            'F': (-0.458, 0.0025),
            'Cl': (-0.584, 0.0016),
            'Br': (-0.538, 0.0022),
            'I': (-0.383, 0.0059),
            'N': (-0.012, 0.0077),
            'Se': (-0.328, 0.0304),
            'Si': (-0.028, 0.0264),
            'Sb': (-0.19, 0.0089),
            'Te': (-0.422, 0.0238),
            'V': (-1.818, 0.0058),
            'Cr': (-2.189, 0.009),
            'Mn': (-1.723, 0.0049),
            'Fe': (-2.227, 0.0082),
            'Co': (-1.518, 0.0055),
            'Ni': (-2.442, 0.0086),
            'W': (-4.988, 0.0176),
            'Mo': (-3.367, 0.0082),
            'H': (-0.181, 0.0013)
        }

        self.no_uncertainties_corrections = {
            'oxide': (-0.733, 0.02),
            'peroxide': (-0.464, 0.0495),
            'superoxide': (-0.153, 0.0549),
            'sulfide': (-0.638, 0.0389),
            'F': (-0.414, 0.0208),
            'Cl': (-0.512, 0.0268),
            'Br': (-0.432, 0.0296),
            'I': (-0.31, 0.0351),
            'N': (-0.35, 0.0621),
            'Se': (-0.358, 0.055),
            'Si': (-0.446, 0.1055),
            'Sb': (-0.278, 0.0781),
            'Te': (-0.446, 0.0653),
            'V': (-1.657, 0.0787),
            'Cr': (-1.889, 0.0704),
            'Mn': (-1.593, 0.0715),
            'Fe': (-2.189, 0.0628),
            'Co': (-1.564, 0.08),
            'Ni': (-2.608, 0.0997),
            'W': (-4.446, 0.1244),
            'Mo': (-3.028, 0.1126),
            'H': (-0.155, 0.0335)
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
