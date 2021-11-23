import json
import os
import unittest

from pymatgen.phonon.bandstructure import (
    PhononBandStructure,
    PhononBandStructureSymmLine,
)
from pymatgen.util.testing import PymatgenTest


class PhononBandStructureSymmLineTest(PymatgenTest):
    def setUp(self):
        with open(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "NaCl_phonon_bandstructure.json"),
            encoding="utf-8",
        ) as f:
            d = json.load(f)
            self.bs = PhononBandStructureSymmLine.from_dict(d)

        with open(
            os.path.join(PymatgenTest.TEST_FILES_DIR, "Si_phonon_bandstructure.json"),
            encoding="utf-8",
        ) as f:
            d = json.load(f)
            self.bs2 = PhononBandStructureSymmLine.from_dict(d)

    def test_basic(self):
        self.assertAlmostEqual(self.bs.bands[1][10], 0.7753555184)
        self.assertAlmostEqual(self.bs.bands[5][100], 5.2548379776)
        self.assertArrayEqual(self.bs.bands.shape, (6, 204))
        self.assertArrayEqual(self.bs.eigendisplacements.shape, (6, 204, 2, 3))
        self.assertArrayAlmostEqual(
            self.bs.eigendisplacements[3][50][0],
            [0.0 + 0.0j, 0.14166569 + 0.04098339j, -0.14166569 - 0.04098339j],
        )
        self.assertTrue(self.bs.has_eigendisplacements, True)

        self.assertArrayEqual(self.bs.min_freq()[0].frac_coords, [0, 0, 0])
        self.assertAlmostEqual(self.bs.min_freq()[1], -0.03700895020)
        self.assertTrue(self.bs.has_imaginary_freq())
        self.assertFalse(self.bs.has_imaginary_freq(tol=0.5))
        self.assertArrayAlmostEqual(self.bs.asr_breaking(), [-0.0370089502, -0.0370089502, -0.0221388897])

        self.assertEqual(self.bs.nb_bands, 6)
        self.assertEqual(self.bs.nb_qpoints, 204)

        self.assertArrayAlmostEqual(self.bs.qpoints[1].frac_coords, [0.01, 0, 0])

    def test_nac(self):
        self.assertTrue(self.bs.has_nac)
        self.assertFalse(self.bs2.has_nac)
        self.assertAlmostEqual(self.bs.get_nac_frequencies_along_dir([1, 1, 0])[3], 4.6084532143)
        self.assertIsNone(self.bs.get_nac_frequencies_along_dir([0, 1, 1]))
        self.assertIsNone(self.bs2.get_nac_frequencies_along_dir([0, 0, 1]))
        self.assertArrayAlmostEqual(
            self.bs.get_nac_eigendisplacements_along_dir([1, 1, 0])[3][1],
            [(0.1063906409128248 + 0j), 0j, 0j],
        )
        self.assertIsNone(self.bs.get_nac_eigendisplacements_along_dir([0, 1, 1]))
        self.assertIsNone(self.bs2.get_nac_eigendisplacements_along_dir([0, 0, 1]))

    def test_branches(self):
        self.assertEqual(self.bs.branches[0]["end_index"], 50)
        self.assertEqual(self.bs.branches[1]["start_index"], 51)
        self.assertEqual(self.bs.branches[2]["name"], "Y-Gamma")
        self.assertAlmostEqual(self.bs.get_branch(10)[0]["name"], "Gamma-X")
        self.assertEqual(len(self.bs.branches), 4)

    def test_dict_methods(self):
        s = self.bs.as_dict()
        self.assertIsNotNone(s)
        self.assertIsNotNone(json.dumps(s))
        s = self.bs2.as_dict()
        self.assertIsNotNone(s)
        self.assertIsNotNone(json.dumps(s))
        s = self.bs2.as_phononwebsite()
        self.assertIsNotNone(s)
        self.assertIsNotNone(json.dumps(s))
        self.assertMSONable(self.bs)
        self.assertMSONable(self.bs2)

    def test_write_methods(self):
        self.bs2.write_phononwebsite("test.json")

    def tearDown(self):
        if os.path.isfile("test.json"):
            os.remove("test.json")


if __name__ == "__main__":
    unittest.main()
