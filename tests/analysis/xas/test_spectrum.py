from __future__ import annotations

import json

import numpy as np
import pytest
from monty.json import MontyDecoder
from numpy.testing import assert_array_equal
from pytest import approx

from pymatgen.analysis.xas.spectrum import XAS, site_weighted_spectrum
from pymatgen.core import Element
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

test_dir = f"{TEST_FILES_DIR}/spectrum_test"

with open(f"{test_dir}/LiCoO2_k_xanes.json") as fp:
    k_xanes_dict = json.load(fp, cls=MontyDecoder)
with open(f"{test_dir}/LiCoO2_k_exafs.json") as fp:
    k_exafs_dict = json.load(fp, cls=MontyDecoder)
with open(f"{test_dir}/ZnO_l2_xanes.json") as fp:
    l2_xanes_dict = json.load(fp, cls=MontyDecoder)
with open(f"{test_dir}/ZnO_l3_xanes.json") as fp:
    l3_xanes_dict = json.load(fp, cls=MontyDecoder)
with open(f"{test_dir}/site1_k_xanes.json") as fp:
    site1_xanes_dict = json.load(fp, cls=MontyDecoder)
with open(f"{test_dir}/site2_k_xanes.json") as fp:
    site2_xanes_dict = json.load(fp, cls=MontyDecoder)


class TestXAS(PymatgenTest):
    def setUp(self):
        self.k_xanes = XAS.from_dict(k_xanes_dict)
        self.k_exafs = XAS.from_dict(k_exafs_dict)
        self.l2_xanes = XAS.from_dict(l2_xanes_dict)
        self.l3_xanes = XAS.from_dict(l3_xanes_dict)
        self.site1_xanes = XAS.from_dict(site1_xanes_dict)
        self.site2_xanes = XAS.from_dict(site2_xanes_dict)

    def test_e0(self):
        assert approx(self.k_xanes.e0) == 7728.565

    def test_k(self):
        assert len(self.k_xanes.x) == len(self.k_xanes.k)
        assert self.k_xanes.e0 == approx(self.k_xanes.x[self.k_xanes.k.index(0)])

    def test_normalization(self):
        self.k_xanes.normalize(mode="sum")
        assert approx(np.sum(self.k_xanes.y)) == 1.0

    def test_add_mul(self):
        scaled_spect = self.k_xanes + self.k_xanes
        scaled_spect2 = self.k_xanes * 3
        assert np.allclose(scaled_spect.y, 2 * self.k_xanes.y)
        assert np.allclose(scaled_spect2.y, 3 * self.k_xanes.y)
        assert approx(self.k_xanes.get_interpolated_value(7720.422), abs=1e-3) == 0.274302

    def test_to_from_dict(self):
        s = XAS.from_dict(self.k_xanes.as_dict())
        assert np.allclose(s.y, self.k_xanes.y)

    def test_attributes(self):
        assert_array_equal(self.k_xanes.energy, self.k_xanes.x)
        assert_array_equal(self.k_xanes.intensity, self.k_xanes.y)

    def test_str(self):
        assert str(self.k_xanes) is not None

    def test_validate(self):
        y_zeros = np.zeros(len(self.k_xanes.x))
        with pytest.raises(ValueError, match="Double check the intensities. Most of them are non-positive"):
            XAS(
                self.k_xanes.x,
                y_zeros,
                self.k_xanes.structure,
                self.k_xanes.absorbing_element,
            )

    def test_stitch_xafs(self):
        with pytest.raises(ValueError, match="Invalid mode. Only XAFS and L23 are supported"):
            XAS.stitch(self.k_xanes, self.k_exafs, mode="invalid")
        xafs = XAS.stitch(self.k_xanes, self.k_exafs, mode="XAFS")
        assert isinstance(xafs, XAS)
        assert xafs.spectrum_type == "XAFS"
        assert len(xafs.x) == 500
        assert min(xafs.x) == approx(min(self.k_xanes.x), abs=1e-2)
        assert max(xafs.y) == approx(max(self.k_xanes.y), abs=1e-2)
        assert xafs.x[np.argmax(np.gradient(xafs.y) / np.gradient(xafs.x))] == approx(self.k_xanes.e0, abs=1e-2)
        with pytest.raises(ValueError, match="The input structures for spectra mismatch"):
            XAS.stitch(self.k_xanes, self.l2_xanes, mode="XAFS")
        self.k_xanes.x = np.zeros(100)
        with pytest.raises(ValueError, match="Energy overlap between XANES and EXAFS is needed for stitching"):
            XAS.stitch(self.k_xanes, self.k_exafs)
        self.k_xanes.absorbing_element = Element("Pt")
        with pytest.raises(ValueError, match="The absorbing elements for spectra are different"):
            XAS.stitch(self.k_xanes, self.k_exafs, mode="XAFS")

    def test_stitch_l23(self):
        self.l2_xanes.y[0] = 0.1
        with pytest.warns(UserWarning, match="jump") as warns:
            XAS.stitch(self.l2_xanes, self.l3_xanes, 100, mode="L23")
            assert len(warns) == 1
        self.l2_xanes = XAS.from_dict(l2_xanes_dict)
        l23 = XAS.stitch(self.l2_xanes, self.l3_xanes, 100, mode="L23")
        assert isinstance(l23, XAS)
        assert l23.edge == "L23"
        assert min(l23.x) == approx(min(self.l3_xanes.x), abs=1e-3)
        assert max(l23.x) == approx(max(self.l3_xanes.x), abs=1e-3)
        assert np.greater_equal(l23.y, self.l2_xanes.y).all()
        assert len(l23.x) == 100
        self.l2_xanes.spectrum_type = "EXAFS"
        with pytest.raises(ValueError, match="Only XANES spectrum can be stitched in L23 mode"):
            XAS.stitch(self.l2_xanes, self.l3_xanes, mode="L23")
        self.l2_xanes.absorbing_element = Element("Pt")
        with pytest.raises(ValueError, match="The absorbing elements for spectra are different"):
            XAS.stitch(self.l2_xanes, self.l3_xanes, mode="L23")
        with pytest.raises(ValueError, match="The input structures for spectra mismatch"):
            XAS.stitch(self.k_xanes, self.l3_xanes, mode="L23")

    def test_site_weighted_spectrum(self):
        weighted_spectrum = site_weighted_spectrum([self.site1_xanes, self.site2_xanes])
        assert isinstance(weighted_spectrum, XAS)
        assert len(weighted_spectrum.x), 500
        # The site multiplicities for site1 and site2 are 4 and 2, respectively.
        assert weighted_spectrum.y[0] == approx((4 * self.site1_xanes.y[0] + 2 * self.site2_xanes.y[0]) / 6, abs=1e-2)
        assert min(weighted_spectrum.x) == max(min(self.site1_xanes.x), min(self.site2_xanes.x))
        self.site2_xanes.absorbing_index = self.site1_xanes.absorbing_index
        with pytest.raises(ValueError, match="Need at least two site-wise spectra to perform site-weighting"):
            site_weighted_spectrum([self.site1_xanes, self.site2_xanes])
