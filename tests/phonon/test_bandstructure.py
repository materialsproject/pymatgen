from __future__ import annotations

import json

from numpy.testing import assert_allclose, assert_array_equal
from pytest import approx

from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestPhononBandStructureSymmLine(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_FILES_DIR}/NaCl_phonon_bandstructure.json") as file:
            dct = json.load(file)
        self.bs = PhononBandStructureSymmLine.from_dict(dct)

        with open(f"{TEST_FILES_DIR}/Si_phonon_bandstructure.json") as file:
            dct = json.load(file)
        self.bs2 = PhononBandStructureSymmLine.from_dict(dct)

    def test_repr(self):
        assert repr(self.bs) == "PhononBandStructureSymmLine(bands=(6, 204), labels=['Gamma', 'X', 'Y', 'Z'])"
        assert repr(self.bs2) == (
            r"PhononBandStructureSymmLine(bands=(6, 130), labels=['$\\Gamma$', 'X', 'W', 'K', 'L', 'U'])"
        )

    def test_basic(self):
        assert self.bs.bands[1][10] == approx(0.7753555184)
        assert self.bs.bands[5][100] == approx(5.2548379776)
        assert_array_equal(self.bs.bands.shape, (6, 204))
        assert_array_equal(self.bs.eigendisplacements.shape, (6, 204, 2, 3))
        assert_allclose(
            self.bs.eigendisplacements[3][50][0],
            [0.0 + 0.0j, 0.14166569 + 0.04098339j, -0.14166569 - 0.04098339j],
        )
        assert self.bs.has_eigendisplacements, True

        assert list(self.bs.min_freq()[0].frac_coords) == [0, 0, 0]
        assert self.bs.min_freq()[1] == approx(-0.03700895020)
        assert_allclose(self.bs.asr_breaking(), [-0.0370089502, -0.0370089502, -0.0221388897])

        assert self.bs.nb_bands == 6
        assert self.bs.nb_qpoints == 204

        assert_allclose(self.bs.qpoints[1].frac_coords, [0.01, 0, 0])

    def test_has_imaginary_freq(self):
        for tol in (0, 0.01, 0.02, 0.03, 0.04, 0.05):
            assert self.bs.has_imaginary_freq(tol=tol) == (tol < 0.04)

        for tol in (0, 0.01, 0.02, 0.03, 0.04, 0.05):
            assert self.bs2.has_imaginary_freq(tol=tol) == (tol < 0.01)

        # test Gamma point imaginary frequency detection
        for tol in (0, 0.01, 0.02, 0.03, 0.04, 0.05):
            assert self.bs.has_imaginary_gamma_freq(tol=tol) == (tol < 0.04)

        for tol in (0, 0.01, 0.02, 0.03, 0.04, 0.05):
            assert self.bs2.has_imaginary_gamma_freq(tol=tol) == (tol < 0.01)

    def test_nac(self):
        assert self.bs.has_nac
        assert not self.bs2.has_nac
        assert self.bs.get_nac_frequencies_along_dir([1, 1, 0])[3] == approx(4.6084532143)
        assert self.bs.get_nac_frequencies_along_dir([0, 1, 1]) is None
        assert self.bs2.get_nac_frequencies_along_dir([0, 0, 1]) is None
        assert_allclose(
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
        dct = self.bs.as_dict()
        assert dct is not None
        assert json.dumps(dct) is not None
        dct = self.bs2.as_dict()
        assert dct is not None
        assert json.dumps(dct) is not None
        dct = self.bs2.as_phononwebsite()
        assert dct is not None
        assert json.dumps(dct) is not None
        self.assert_msonable(self.bs)
        self.assert_msonable(self.bs2)

    def test_write_methods(self):
        self.bs2.write_phononwebsite(f"{self.tmp_path}/test.json")
