from __future__ import annotations

import json
import os

import numpy as np
from numpy.testing import assert_array_equal
from pytest import approx

from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestPhononBandStructureSymmLine(PymatgenTest):
    def setUp(self):
        with open(
            f"{TEST_FILES_DIR}/NaCl_phonon_bandstructure.json",
            encoding="utf-8",
        ) as file:
            dct = json.load(file)
            self.bs = PhononBandStructureSymmLine.from_dict(dct)

        with open(
            f"{TEST_FILES_DIR}/Si_phonon_bandstructure.json",
            encoding="utf-8",
        ) as file:
            dct = json.load(file)
            self.bs2 = PhononBandStructureSymmLine.from_dict(dct)

    def test_basic(self):
        assert self.bs.bands[1][10] == approx(0.7753555184)
        assert self.bs.bands[5][100] == approx(5.2548379776)
        assert_array_equal(self.bs.bands.shape, (6, 204))
        assert_array_equal(self.bs.eigendisplacements.shape, (6, 204, 2, 3))
        assert np.allclose(
            self.bs.eigendisplacements[3][50][0],
            [0.0 + 0.0j, 0.14166569 + 0.04098339j, -0.14166569 - 0.04098339j],
        )
        assert self.bs.has_eigendisplacements, True

        assert_array_equal(self.bs.min_freq()[0].frac_coords, [0, 0, 0])
        assert self.bs.min_freq()[1] == approx(-0.03700895020)
        assert self.bs.has_imaginary_freq()
        assert not self.bs.has_imaginary_freq(tol=0.5)
        assert np.allclose(self.bs.asr_breaking(), [-0.0370089502, -0.0370089502, -0.0221388897])

        assert self.bs.nb_bands == 6
        assert self.bs.nb_qpoints == 204

        assert np.allclose(self.bs.qpoints[1].frac_coords, [0.01, 0, 0])

    def test_nac(self):
        assert self.bs.has_nac
        assert not self.bs2.has_nac
        assert self.bs.get_nac_frequencies_along_dir([1, 1, 0])[3] == approx(4.6084532143)
        assert self.bs.get_nac_frequencies_along_dir([0, 1, 1]) is None
        assert self.bs2.get_nac_frequencies_along_dir([0, 0, 1]) is None
        assert np.allclose(
            self.bs.get_nac_eigendisplacements_along_dir([1, 1, 0])[3][1],
            [(0.1063906409128248 + 0j), 0j, 0j],
        )
        assert self.bs.get_nac_eigendisplacements_along_dir([0, 1, 1]) is None
        assert self.bs2.get_nac_eigendisplacements_along_dir([0, 0, 1]) is None

    def test_branches(self):
        assert self.bs.branches[0]["end_index"] == 50
        assert self.bs.branches[1]["start_index"] == 51
        assert self.bs.branches[2]["name"] == "Y-Gamma"
        assert self.bs.get_branch(10)[0]["name"] == "Gamma-X"
        assert len(self.bs.branches) == 4

    def test_dict_methods(self):
        s = self.bs.as_dict()
        assert s is not None
        assert json.dumps(s) is not None
        s = self.bs2.as_dict()
        assert s is not None
        assert json.dumps(s) is not None
        s = self.bs2.as_phononwebsite()
        assert s is not None
        assert json.dumps(s) is not None
        self.assert_msonable(self.bs)
        self.assert_msonable(self.bs2)

    def test_write_methods(self):
        self.bs2.write_phononwebsite("test.json")

    def tearDown(self):
        if os.path.isfile("test.json"):
            os.remove("test.json")
