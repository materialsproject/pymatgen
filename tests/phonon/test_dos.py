from __future__ import annotations

import json
import re

import numpy as np
import pytest
from pytest import approx

from pymatgen.core import Element
from pymatgen.phonon.dos import CompletePhononDos, PhononDos
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

TEST_DIR = f"{TEST_FILES_DIR}/phonon/dos"


class TestPhononDos(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_DIR}/NaCl_ph_dos.json") as file:
            self.dos = PhononDos.from_dict(json.load(file))
        with open(f"{TEST_DIR}/NaCl_complete_ph_dos.json") as file:
            self.structure = CompletePhononDos.from_dict(json.load(file)).structure

    def test_repr(self):
        assert repr(self.dos) == "PhononDos(frequencies=(201,), densities=(201,), n_positive_freqs=183)"

    def test_str(self):
        assert re.match(
            r"#Frequency\s+Density\s+\n-0.66954\s+0.00000\n-0.63158\s+0.00000\n-0.59363\s+0.00000", str(self.dos)
        )

    def test_properties(self):
        assert self.dos.densities[15] == approx(0.0001665998)
        assert self.dos.frequencies[20] == approx(0.0894965119)
        assert self.dos.get_interpolated_value(3.0) == approx(1.2915532670115628)
        assert len(self.dos.frequencies) == 201
        assert len(self.dos.densities) == 201

    def test_get_smeared_densities(self):
        smeared = self.dos.get_smeared_densities(0.01)
        assert smeared[20] == approx(0.00084171007635058825)
        dens = self.dos.densities
        assert sum(dens) == approx(sum(smeared))

        # test 0 smearing returns original DOS
        assert self.dos.get_smeared_densities(0) is self.dos.densities

    def test_dict_methods(self):
        json_str = json.dumps(self.dos.as_dict())
        assert json_str.startswith('{"@module": "pymatgen.phonon.dos", "@class": "PhononDos", "frequencies":')
        self.assert_msonable(self.dos)

    def test_thermodynamic_functions(self):
        assert self.dos.cv(300, structure=self.structure) == approx(48.049366665412485, abs=1e-4)
        assert self.dos.internal_energy(300, structure=self.structure) == approx(15527.596956593827, abs=1e-4)
        assert self.dos.helmholtz_free_energy(300, structure=self.structure) == approx(-6998.034212172695, abs=1e-4)
        assert self.dos.entropy(300, structure=self.structure) == approx(75.08543723748751, abs=1e-4)
        assert self.dos.zero_point_energy(structure=self.structure) == approx(4847.462485708741, abs=1e-4)

    def test_add(self):
        dos_2x = self.dos + self.dos
        assert dos_2x.frequencies == approx(self.dos.frequencies)
        assert dos_2x.densities == approx(2 * self.dos.densities)

        dos_3x = self.dos + dos_2x
        assert dos_3x.frequencies == approx(self.dos.frequencies)
        assert dos_3x.densities == approx(3 * self.dos.densities)

        # test commutativity
        assert dos_2x + 42 == 42 + dos_2x

    def test_sub(self):
        dos_0 = self.dos - self.dos
        assert dos_0.frequencies == approx(self.dos.frequencies)
        assert dos_0.densities == approx(self.dos.densities * 0)

        dos_1 = self.dos - dos_0
        assert dos_1.frequencies == approx(self.dos.frequencies)
        assert dos_1.densities == approx(self.dos.densities)

    def test_mul(self):
        dos_2x = self.dos * 2
        assert dos_2x.frequencies == approx(self.dos.frequencies)
        assert dos_2x.densities == approx(2 * self.dos.densities)

        # test commutativity
        assert dos_2x * 1.234 == 1.234 * dos_2x

    def test_eq(self):
        assert self.dos == self.dos
        assert self.dos != 42
        assert self.dos != 2 * self.dos
        assert 2 * self.dos == self.dos + self.dos

    def test_mae(self):
        assert self.dos.mae(self.dos) == 0
        assert self.dos.mae(self.dos + 1) == 1
        assert self.dos.mae(self.dos - 1) == 1
        assert self.dos.mae(2 * self.dos) == pytest.approx(0.786546967)
        assert (2 * self.dos).mae(self.dos) == pytest.approx(0.786546967)

        # test two_sided=False after shifting DOS freqs so MAE requires interpolation
        dos2 = PhononDos(self.dos.frequencies + 0.01, self.dos.densities)
        assert self.dos.mae(dos2 + 1, two_sided=False) == pytest.approx(0.999999999)
        assert self.dos.mae(dos2 - 1, two_sided=False) == pytest.approx(1.00000000000031)
        assert self.dos.mae(2 * dos2, two_sided=False) == pytest.approx(0.786546967)

    def test_r2_score(self):
        assert self.dos.r2_score(self.dos) == 1
        assert self.dos.r2_score(self.dos + 1) == pytest.approx(-0.45647319)
        assert self.dos.r2_score(self.dos - 1) == pytest.approx(-0.45647319)
        assert self.dos.r2_score(2 * self.dos) == pytest.approx(-0.901056070)

        # check that r2_score is 0 for DOS with same mean as self.dos
        densities = self.dos.densities
        mean_dos = PhononDos(self.dos.frequencies, np.full_like(densities, densities.mean()))
        assert self.dos.r2_score(mean_dos) == pytest.approx(0)
        # moving away from the mean should decrease r2_score
        assert self.dos.r2_score(-mean_dos) == pytest.approx(-3.604224283)

    def test_get_last_peak(self):
        peak_freq = self.dos.get_last_peak()
        assert peak_freq == approx(5.9909763191)
        assert self.dos.get_interpolated_value(peak_freq) == approx(1.1700016497)

        # try different threshold
        peak_freq = self.dos.get_last_peak(threshold=0.5)
        assert peak_freq == approx(4.9662820761)

    def test_get_dos_fp(self):
        # normalize is True
        dos_fp = self.dos.get_dos_fp(min_f=-1, max_f=5, n_bins=56, normalize=True)
        bin_width = np.diff(dos_fp.frequencies)[0][0]
        assert max(dos_fp.frequencies[0]) <= 5
        assert min(dos_fp.frequencies[0]) >= -1
        assert len(dos_fp.frequencies[0]) == 56
        assert sum(dos_fp.densities * bin_width) == approx(1)
        # normalize is False
        dos_fp2 = self.dos.get_dos_fp(min_f=-1, max_f=5, n_bins=56, normalize=False)
        bin_width2 = np.diff(dos_fp2.frequencies)[0][0]
        assert sum(dos_fp2.densities * bin_width2) == approx(13.722295798242834)
        assert dos_fp2.bin_width == approx(bin_width2)
        # binning is False
        dos_fp = self.dos.get_dos_fp(min_f=None, max_f=None, n_bins=56, normalize=True, binning=False)
        assert dos_fp.n_bins == len(self.dos.frequencies)

    def test_get_dos_fp_similarity(self):
        # Tanimoto
        dos_fp = self.dos.get_dos_fp(min_f=-1, max_f=6, n_bins=56, normalize=True)
        dos_fp2 = self.dos.get_dos_fp(min_f=-1, max_f=6, n_bins=56, normalize=False)
        similarity_index = self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric="tanimoto")
        assert similarity_index == approx(0.0553088193)

        dos_fp = self.dos.get_dos_fp(min_f=-1, max_f=6, n_bins=56, normalize=True)
        dos_fp2 = self.dos.get_dos_fp(min_f=-1, max_f=6, n_bins=56, normalize=True)
        similarity_index = self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric="tanimoto")
        assert similarity_index == approx(1)

        # Wasserstein
        dos_fp = self.dos.get_dos_fp(min_f=-1, max_f=6, n_bins=56, normalize=True)
        dos_fp2 = self.dos.get_dos_fp(min_f=-1, max_f=6, n_bins=56, normalize=True)
        similarity_index = self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric="wasserstein")
        assert similarity_index == approx(0)

    def test_dos_fp_exceptions(self):
        dos_fp = self.dos.get_dos_fp(min_f=-1, max_f=5, n_bins=56, normalize=True)
        dos_fp2 = self.dos.get_dos_fp(min_f=-1, max_f=5, n_bins=56, normalize=True)
        # test exceptions
        with pytest.raises(
            ValueError,
            match="Cannot compute similarity index. When normalize=True, then please set metric=cosine-sim",
        ):
            self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric="tanimoto", normalize=True)

        valid_metrics = ("tanimoto", "wasserstein", "cosine-sim")
        metric = "Dot"
        with pytest.raises(ValueError, match=re.escape(f"Invalid {metric=}, choose from {valid_metrics}.")):
            self.dos.get_dos_fp_similarity(dos_fp, dos_fp2, col=1, metric=metric, normalize=False)


class TestCompletePhononDos(PymatgenTest):
    def setUp(self):
        with open(f"{TEST_DIR}/NaCl_complete_ph_dos.json") as file:
            self.cdos = CompletePhononDos.from_dict(json.load(file))

    def test_properties(self):
        site_Na = self.cdos.structure[0]
        site_Cl = self.cdos.structure[1]

        assert len(self.cdos.frequencies) == 201
        assert self.cdos.pdos[site_Na][30] == approx(0.008058208)
        assert self.cdos.get_site_dos(site_Na).densities[30] == approx(0.008058208)
        assert self.cdos.pdos[site_Cl][30] == approx(0.0119040783)

        assert Element.Na in self.cdos.get_element_dos()
        assert Element.Cl in self.cdos.get_element_dos()

        sum_dos = self.cdos.get_element_dos()[Element.Na] + self.cdos.get_element_dos()[Element.Cl]
        assert sum_dos.frequencies == approx(self.cdos.frequencies)
        assert sum_dos.densities == approx(self.cdos.densities)

    def test_dict_methods(self):
        json_str = json.dumps(self.cdos.as_dict())
        assert json_str is not None
        self.assert_msonable(self.cdos)

    def test_str(self):
        assert str(self.cdos).startswith(
            "Complete phonon DOS for Full Formula (Na1 Cl1)\nReduced Formula: NaCl\n"
            "abc   :   4.023651   4.023651   4.023651\nangles"
        )
