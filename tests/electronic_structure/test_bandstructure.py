from __future__ import annotations

import json
import unittest
import warnings

import numpy as np
import pytest
from monty.serialization import loadfn
from pytest import approx

from pymatgen.core.lattice import Lattice
from pymatgen.electronic_structure.bandstructure import (
    BandStructureSymmLine,
    Kpoint,
    LobsterBandStructureSymmLine,
    get_reconstructed_band_structure,
)
from pymatgen.electronic_structure.core import Orbital, Spin
from pymatgen.electronic_structure.plotter import BSPlotterProjected
from pymatgen.io.vasp import BSVasprun
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest


class TestKpoint(unittest.TestCase):
    def setUp(self):
        self.lattice = Lattice.cubic(10.0)
        self.kpoint = Kpoint([0.1, 0.4, -0.5], self.lattice, label="X")

    def test_properties(self):
        assert list(self.kpoint.frac_coords) == [0.1, 0.4, -0.5]
        assert self.kpoint.a == 0.1
        assert self.kpoint.b == 0.4
        assert self.kpoint.c == -0.5
        assert self.lattice == Lattice.cubic(10.0)
        assert list(self.kpoint.cart_coords) == [1.0, 4.0, -5.0]
        assert self.kpoint.label == "X"

    def test_as_dict(self):
        assert isinstance(self.kpoint.as_dict()["fcoords"], list)
        assert isinstance(self.kpoint.as_dict()["ccoords"], list)
        assert not isinstance(self.kpoint.as_dict()["fcoords"][0], np.float64)
        assert not isinstance(self.kpoint.as_dict()["ccoords"][0], np.float64)
        assert self.kpoint.as_dict()["fcoords"] == [0.1, 0.4, -0.5]
        assert self.kpoint.as_dict()["ccoords"] == [1.0, 4.0, -5.0]

    def test_from_dict(self):
        d = self.kpoint.as_dict()

        kpoint = Kpoint.from_dict(d)

        assert list(kpoint.frac_coords) == [0.1, 0.4, -0.5]
        assert kpoint.a == 0.1
        assert kpoint.b == 0.4
        assert kpoint.c == -0.5
        assert kpoint.lattice == Lattice.cubic(10.0)
        assert list(kpoint.cart_coords) == [1.0, 4.0, -5.0]
        assert kpoint.label == "X"


class TestBandStructureSymmLine(PymatgenTest):
    def setUp(self):
        self.bs: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/Cu2O_361_bandstructure.json")
        self.bs2: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/CaO_2605_bandstructure.json")
        self.bs_spin: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/NiO_19009_bandstructure.json")
        self.bs_cbm0: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/InN_22205_bandstructure.json")
        self.bs_cu: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/Cu_30_bandstructure.json")
        self.bs_diff_spins: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/VBr2_971787_bandstructure.json")
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_basic(self):
        assert np.allclose(self.bs.projections[Spin.up][10][12][0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        assert np.allclose(
            self.bs.projections[Spin.up][25][0][Orbital.dyz.value],
            [0.0, 0.0, 0.0011, 0.0219, 0.0219, 0.069],
        )
        assert self.bs.get_projection_on_elements()[Spin.up][25][10]["O"] == approx(0.0328)
        assert self.bs.get_projection_on_elements()[Spin.up][22][25]["Cu"] == approx(0.8327)
        proj = self.bs.get_projections_on_elements_and_orbitals({"Cu": ["s", "d"]})
        assert proj[Spin.up][25][0]["Cu"]["s"] == approx(0.0027)
        assert proj[Spin.up][25][0]["Cu"]["d"] == approx(0.8495999999999999)

        assert self.bs2.nb_bands == 16
        assert self.bs2.bands[Spin.up][5][10] == approx(0.5608)
        assert self.bs2.bands[Spin.up][5][10] == approx(0.5608)
        assert self.bs2.branches[5]["name"] == "L-U"
        assert self.bs2.branches[5]["start_index"] == 80
        assert self.bs2.branches[5]["end_index"] == 95
        assert self.bs2.distance[70] == approx(4.2335127528765737)

        assert self.bs_spin.nb_bands == 27
        assert self.bs_spin.bands[Spin.up][5][10] == approx(0.262)
        assert self.bs_spin.bands[Spin.down][5][10] == approx(1.6156)

    def test_properties(self):
        self.one_kpoint = self.bs2.kpoints[31]
        assert list(self.one_kpoint.frac_coords) == [0.5, 0.25, 0.75]
        assert self.one_kpoint.cart_coords == approx([0.64918757, 1.29837513, 0.0])
        assert self.one_kpoint.label == "W"

        assert self.bs2.efermi == approx(2.6211967), "wrong fermi energy"

    def test_get_branch(self):
        assert self.bs2.get_branch(110)[0]["name"] == "U-W"

    def test_get_direct_band_gap_dict(self):
        direct_dict = self.bs_diff_spins.get_direct_band_gap_dict()
        assert direct_dict[Spin.down]["value"] == 4.5365

        for bs in [self.bs2, self.bs_spin]:
            dg_dict = bs.get_direct_band_gap_dict()
            for spin, v in bs.bands.items():
                kpt = dg_dict[spin]["kpoint_index"]
                vb, cb = dg_dict[spin]["band_indices"]
                gap = v[cb][kpt] - v[vb][kpt]
                assert gap == dg_dict[spin]["value"]
        with pytest.raises(ValueError, match="get_direct_band_gap_dict should only be used with non-metals"):
            self.bs_cu.get_direct_band_gap_dict()

    def test_get_direct_band_gap(self):
        assert self.bs2.get_direct_band_gap() == approx(4.0125999999999999)
        assert self.bs_diff_spins.get_direct_band_gap() > 0
        assert self.bs_cu.get_direct_band_gap() == 0

    def test_is_metal(self):
        assert not self.bs2.is_metal(), "wrong metal assignment"
        assert not self.bs_spin.is_metal(), "wrong metal assignment"
        assert self.bs_cu.is_metal(), "wrong metal assignment"

    def test_get_cbm(self):
        cbm = self.bs2.get_cbm()
        assert cbm["energy"] == approx(5.8709), "wrong CBM energy"
        assert cbm["band_index"][Spin.up][0] == 8, "wrong CBM band index"
        assert cbm["kpoint_index"][0] == 15, "wrong CBM kpoint index"
        assert cbm["kpoint"].frac_coords[0] == 0.5, "wrong CBM kpoint frac coords"
        assert cbm["kpoint"].frac_coords[1] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm["kpoint"].frac_coords[2] == 0.5, "wrong CBM kpoint frac coords"
        assert cbm["kpoint"].label == "X", "wrong CBM kpoint label"
        cbm_spin = self.bs_spin.get_cbm()
        assert cbm_spin["energy"] == approx(8.0458), "wrong CBM energy"
        assert cbm_spin["band_index"][Spin.up][0] == 12, "wrong CBM band index"
        assert len(cbm_spin["band_index"][Spin.down]) == 0, "wrong CBM band index"
        assert cbm_spin["kpoint_index"][0] == 0, "wrong CBM kpoint index"
        assert cbm_spin["kpoint"].frac_coords[0] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm_spin["kpoint"].frac_coords[1] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm_spin["kpoint"].frac_coords[2] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm_spin["kpoint"].label == "\\Gamma", "wrong CBM kpoint label"

    def test_get_vbm(self):
        vbm = self.bs2.get_vbm()
        assert vbm["energy"] == approx(2.2361), "wrong VBM energy"
        assert len(vbm["band_index"][Spin.up]) == 3, "wrong VBM number of bands"
        assert vbm["band_index"][Spin.up][0] == 5, "wrong VBM band index"
        assert vbm["kpoint_index"][0] == 0, "wrong VBM kpoint index"
        assert vbm["kpoint"].frac_coords[0] == 0.0, "wrong VBM kpoint frac coords"
        assert vbm["kpoint"].frac_coords[1] == 0.0, "wrong VBM kpoint frac coords"
        assert vbm["kpoint"].frac_coords[2] == 0.0, "wrong VBM kpoint frac coords"
        assert vbm["kpoint"].label == "\\Gamma", "wrong VBM kpoint label"
        vbm_spin = self.bs_spin.get_vbm()
        assert vbm_spin["energy"] == approx(5.731), "wrong VBM energy"
        assert len(vbm_spin["band_index"][Spin.up]) == 2, "wrong VBM number of bands"
        assert len(vbm_spin["band_index"][Spin.down]) == 0, "wrong VBM number of bands"
        assert vbm_spin["band_index"][Spin.up][0] == 10, "wrong VBM band index"
        assert vbm_spin["kpoint_index"][0] == 79, "wrong VBM kpoint index"
        assert vbm_spin["kpoint"].frac_coords[0] == 0.5, "wrong VBM kpoint frac coords"
        assert vbm_spin["kpoint"].frac_coords[1] == 0.5, "wrong VBM kpoint frac coords"
        assert vbm_spin["kpoint"].frac_coords[2] == 0.5, "wrong VBM kpoint frac coords"
        assert vbm_spin["kpoint"].label == "L", "wrong VBM kpoint label"

    def test_get_band_gap(self):
        bg = self.bs2.get_band_gap()
        assert bg["energy"] == approx(3.6348), "wrong gap energy"
        assert bg["transition"] == "\\Gamma-X", "wrong kpoint transition"
        assert not bg["direct"], "wrong nature of the gap"
        bg_spin = self.bs_spin.get_band_gap()
        assert bg_spin["energy"] == approx(2.3148), "wrong gap energy"
        assert bg_spin["transition"] == "L-\\Gamma", "wrong kpoint transition"
        assert not bg_spin["direct"], "wrong nature of the gap"
        bg_cbm0 = self.bs_cbm0.get_band_gap()
        assert bg_cbm0["energy"] == approx(0, abs=1e-3), "wrong gap energy"

    def test_get_sym_eq_kpoints_and_degeneracy(self):
        bs = self.bs2
        cbm_k = bs.get_cbm()["kpoint"].frac_coords
        vbm_k = bs.get_vbm()["kpoint"].frac_coords
        assert bs.get_kpoint_degeneracy(cbm_k) is None
        bs.structure: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/CaO_2605_structure.json")
        assert bs.get_kpoint_degeneracy(cbm_k) == 3
        assert bs.get_kpoint_degeneracy(vbm_k) == 1
        cbm_eqs = bs.get_sym_eq_kpoints(cbm_k)
        assert [0.5, 0.0, 0.5] in cbm_eqs
        assert [0.0, 0.5, 0.5] in cbm_eqs
        assert [0.5, 0.5, 0.0] in cbm_eqs
        vbm_eqs = bs.get_sym_eq_kpoints(vbm_k)
        assert [0.0, 0.0, 0.0] in vbm_eqs

    def test_as_dict(self):
        expected_keys = {
            "@module",
            "@class",
            "lattice_rec",
            "efermi",
            "kpoints",
            "bands",
            "is_metal",
            "vbm",
            "cbm",
            "band_gap",
            "labels_dict",
            "is_spin_polarized",
            "projections",
            # "structure",  # not always present
            "branches",
        }
        d1 = self.bs.as_dict()
        assert set(d1) >= expected_keys, f"{expected_keys - set(d1)=}"
        d2 = self.bs2.as_dict()
        assert set(d2) >= expected_keys, f"{expected_keys - set(d2)=}"
        d3 = self.bs_spin.as_dict()
        assert set(d3) >= expected_keys, f"{expected_keys - set(d3)=}"

    def test_old_format_load(self):
        with open(f"{TEST_FILES_DIR}/bs_ZnS_old.json") as f:
            d = json.load(f)
            bs_old = BandStructureSymmLine.from_dict(d)
            assert bs_old.get_projection_on_elements()[Spin.up][0][0]["Zn"] == 0.0971

    def test_apply_scissor_insulator(self):
        # test applying a scissor operator to a metal
        for scissor in (1, 3):
            bs_scissored = self.bs.apply_scissor(scissor)
            assert not bs_scissored.is_metal()
            assert bs_scissored.nb_bands == 48
            assert bs_scissored.efermi == approx(3.75640309 + scissor)
            orig_efermi = self.bs_spin.efermi
            assert bs_scissored.efermi != approx(orig_efermi)

    def test_apply_scissor_spin_polarized(self):
        # test applying a scissor operator to a spin-polarized system
        bs_scissored = self.bs_spin.apply_scissor(1.0)
        assert bs_scissored.is_metal()
        assert bs_scissored.nb_bands == 27
        assert {*bs_scissored.bands} == {Spin.up, Spin.down}
        assert bs_scissored.efermi == approx(4.64005999)
        orig_efermi = self.bs_spin.efermi
        assert bs_scissored.efermi != approx(orig_efermi)


class TestReconstructBandStructure(PymatgenTest):
    def setUp(self):
        self.bs_cu: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/Cu_30_bandstructure.json")
        self.bs_cu2: BandStructureSymmLine = loadfn(f"{TEST_FILES_DIR}/Cu_30_bandstructure.json")
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def test_reconstruct_band_structure(self):
        bs = get_reconstructed_band_structure([self.bs_cu, self.bs_cu2])
        assert bs.bands[Spin.up].shape == (20, 700), "wrong number of bands or kpoints"

    def test_vasprun_bs(self):
        bsv = BSVasprun(
            f"{TEST_FILES_DIR}/vasprun.xml",
            parse_projected_eigen=True,
            parse_potcar_file=True,
        )
        bs = bsv.get_band_structure(kpoints_filename=f"{TEST_FILES_DIR}/KPOINTS.band", line_mode=True)
        bs.get_projection_on_elements()


class TestLobsterBandStructureSymmLine(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")
        with open(
            f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p/lobster_band_structure_spin.json",
        ) as f:
            bs_spin_dict = json.load(f)
        self.bs_spin = LobsterBandStructureSymmLine.from_dict(bs_spin_dict)

        with open(
            f"{TEST_FILES_DIR}/cohp/Fatband_SiO2/Test_p/lobster_band_structure.json",
        ) as f:
            bs_dict = json.load(f)
        self.bs_p = LobsterBandStructureSymmLine.from_dict(bs_dict)

    def tearDown(self):
        warnings.simplefilter("default")

    def test_basic(self):
        bs_p = self.bs_p
        bs_spin = self.bs_spin
        assert bs_p.structure[0].frac_coords == approx([0.0, 0.47634315, 0.666667])
        assert bs_p.structure[0].species_string == "Si"
        assert bs_p.structure[0].coords == approx([-1.19607309, 2.0716597, 3.67462144])
        assert bs_p.efermi == approx(1.06470288)

        lattice = bs_p.lattice_rec.as_dict()
        assert lattice["matrix"][0] == approx([1.2511575194890285, 0.7223560132915973, 0.0])
        assert lattice["matrix"][1] == approx([0.0, 1.4447123171425553, 0.0])
        assert lattice["matrix"][2] == approx([0.0, 0.0, 1.1399248502312707])
        assert bs_p.kpoints[8].frac_coords == approx([0.09090909, 0.0, 0.0])
        assert bs_p.kpoints[8].cart_coords == approx([0.11374159, 0.06566873, 0.0])
        assert bs_p.kpoints[50].frac_coords == approx([0.46153846, 0.07692308, 0.0])
        assert bs_p.kpoints[50].cart_coords == approx([0.57745732, 0.4445268, 0.0])
        assert bs_p.distance[30] == approx(0.49251552363382556)
        assert bs_p.branches[0]["name"], "\\Gamma-K"
        assert bs_p.get_band_gap()["energy"] == approx(5.6739999999999995)
        assert bs_p.get_projection_on_elements()[Spin.up][0][0]["Si"] == approx(3 * (0.001 + 0.064))
        assert bs_p.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.up][0][0]["Si"]["3p"] == approx(0.003)
        assert bs_p.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.up][0][0]["O"]["2p"] == approx(
            0.002 * 3 + 0.003 * 3
        )
        dict_here = bs_p.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[Spin.up][0][
            0
        ]
        assert dict_here["Si"]["3s"] == approx(0.192)
        assert dict_here["Si"]["3p"] == approx(0.003)
        assert dict_here["O"]["2s"] == approx(0.792)
        assert dict_here["O"]["2p"] == approx(0.015)

        assert bs_spin.get_projection_on_elements()[Spin.up][0][0]["Si"] == approx(3 * (0.001 + 0.064))
        assert bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.up][0][0]["Si"]["3p"] == approx(
            0.003
        )
        assert bs_spin.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.up][0][0]["O"]["2p"] == approx(
            0.002 * 3 + 0.003 * 3
        )

        dict_here = bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[Spin.up][
            0
        ][0]
        assert dict_here["Si"]["3s"] == approx(0.192)
        assert dict_here["Si"]["3p"] == approx(0.003)
        assert dict_here["O"]["2s"] == approx(0.792)
        assert dict_here["O"]["2p"] == approx(0.015)
        assert bs_spin.get_projection_on_elements()[Spin.up][0][0]["Si"] == approx(3 * (0.001 + 0.064))
        assert bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3p"]})[Spin.down][0][0]["Si"]["3p"] == approx(
            0.003
        )
        assert bs_spin.get_projections_on_elements_and_orbitals({"O": ["2p"]})[Spin.down][0][0]["O"]["2p"] == approx(
            0.002 * 3 + 0.003 * 3
        )
        dict_here = bs_spin.get_projections_on_elements_and_orbitals({"Si": ["3s", "3p"], "O": ["2s", "2p"]})[
            Spin.down
        ][0][0]
        assert dict_here["Si"]["3s"] == approx(0.192)
        assert dict_here["Si"]["3p"] == approx(0.003)
        assert dict_here["O"]["2s"] == approx(0.792)
        assert dict_here["O"]["2p"] == approx(0.015)

    def test_proj_bandstructure_plot(self):
        # make sure that it can be plotted!
        BSPlotterProjected(self.bs_spin).get_elt_projected_plots()
        BSPlotterProjected(self.bs_spin).get_projected_plots_dots({"Si": ["3s"]})

    def test_get_branch(self):
        branch = self.bs_p.get_branch(0)[0]
        assert branch["name"] == "\\Gamma-K"
        assert branch["start_index"] == 0
        assert branch["end_index"] == 70
        assert branch["index"] == 0

    def test_get_direct_band_gap_dict(self):
        direct_dict = self.bs_p.get_direct_band_gap_dict()
        assert direct_dict[Spin.up]["value"] == approx(6.005999999999999)
        assert direct_dict[Spin.up]["kpoint_index"] == 0
        assert direct_dict[Spin.up]["band_indices"] == [22, 24]

        direct_dict = self.bs_spin.get_direct_band_gap_dict()
        assert direct_dict[Spin.up]["value"] == approx(6.005999999999999)
        assert direct_dict[Spin.up]["kpoint_index"] == 0
        assert direct_dict[Spin.up]["band_indices"] == [22, 24]
        assert direct_dict[Spin.down]["value"] == approx(6.005999999999999)
        assert direct_dict[Spin.down]["kpoint_index"] == 0
        assert direct_dict[Spin.down]["band_indices"] == [22, 24]

    def test_get_direct_band_gap(self):
        assert self.bs_p.get_direct_band_gap() == approx(6.005999999999999)
        assert self.bs_spin.get_direct_band_gap() == approx(6.005999999999999)

    def test_is_metal(self):
        assert not self.bs_p.is_metal(), "wrong metal assignment"
        assert not self.bs_spin.is_metal(), "wrong metal assignment"

    def test_get_cbm(self):
        cbm = self.bs_p.get_cbm()
        assert cbm["energy"] == approx(6.3037028799999995), "wrong CBM energy"
        assert cbm["band_index"][Spin.up][0] == 24, "wrong CBM band index"
        assert cbm["kpoint_index"][0] == 0, "wrong CBM kpoint index"
        assert cbm["kpoint"].frac_coords[0] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm["kpoint"].frac_coords[1] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm["kpoint"].frac_coords[2] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm["kpoint"].label == "\\Gamma", "wrong CBM kpoint label"
        cbm_spin = self.bs_spin.get_cbm()
        assert cbm_spin["energy"] == approx(6.30370274), "wrong CBM energy"
        assert cbm_spin["band_index"][Spin.up][0] == 24, "wrong CBM band index"
        assert len(cbm_spin["band_index"][Spin.down]) == 1, "wrong CBM band index"
        assert cbm_spin["kpoint_index"][0] == 0, "wrong CBM kpoint index"
        assert cbm_spin["kpoint"].frac_coords[0] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm_spin["kpoint"].frac_coords[1] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm_spin["kpoint"].frac_coords[2] == 0.0, "wrong CBM kpoint frac coords"
        assert cbm_spin["kpoint"].label == "\\Gamma", "wrong CBM kpoint label"

    def test_get_vbm(self):
        vbm = self.bs_p.get_vbm()
        assert vbm["energy"] == approx(0.62970288), "wrong VBM energy"
        assert len(vbm["band_index"][Spin.up]) == 1, "wrong VBM number of bands"
        assert vbm["band_index"][Spin.up][0] == 23, "wrong VBM band index"
        assert vbm["kpoint_index"][0] == 68, "wrong VBM kpoint index"
        assert vbm["kpoint"].frac_coords == approx(
            [0.34615384615385, 0.30769230769231, 0.0]
        ), "wrong VBM kpoint frac coords"
        assert vbm["kpoint"].label is None, "wrong VBM kpoint label"
        vbm_spin = self.bs_spin.get_vbm()
        assert vbm_spin["energy"] == approx(0.6297027399999999), "wrong VBM energy"
        assert len(vbm_spin["band_index"][Spin.up]) == 1, "wrong VBM number of bands"
        assert len(vbm_spin["band_index"][Spin.down]) == 1, "wrong VBM number of bands"
        assert vbm_spin["band_index"][Spin.up][0] == 23, "wrong VBM band index"
        assert vbm_spin["kpoint_index"][0] == 68, "wrong VBM kpoint index"
        assert vbm_spin["kpoint"].frac_coords == approx(
            [0.34615384615385, 0.30769230769231, 0.0]
        ), "wrong VBM kpoint frac coords"
        assert vbm_spin["kpoint"].label is None, "wrong VBM kpoint label"

    def test_get_band_gap(self):
        bg = self.bs_p.get_band_gap()
        assert bg["energy"] == approx(5.6739999999999995), "wrong gap energy"
        assert bg["transition"] == "(0.346,0.308,0.000)-\\Gamma", "wrong kpoint transition"
        assert not bg["direct"], "wrong nature of the gap"
        bg_spin = self.bs_spin.get_band_gap()
        assert bg_spin["energy"] == approx(5.674), "wrong gap energy"
        assert bg_spin["transition"] == "(0.346,0.308,0.000)-\\Gamma", "wrong kpoint transition"
        assert not bg_spin["direct"], "wrong nature of the gap"

    def test_get_sym_eq_kpoints_and_degeneracy(self):
        bs = self.bs_p
        cbm_k = bs.get_cbm()["kpoint"].frac_coords
        vbm_k = bs.get_vbm()["kpoint"].frac_coords
        assert bs.get_kpoint_degeneracy(cbm_k) == 1
        assert bs.get_kpoint_degeneracy(vbm_k) == 3

    def test_as_dict(self):
        dict_str = json.dumps(self.bs_p.as_dict())
        assert dict_str is not None
        dict_str = json.dumps(self.bs_spin.as_dict())
        assert dict_str is not None

    def test_old_format_load(self):
        # this method will use the loading from the old dict
        self.bs_spin.apply_scissor(3.0)
