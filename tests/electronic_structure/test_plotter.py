from __future__ import annotations

import os
from shutil import which

import matplotlib.pyplot as plt
import numpy as np
import orjson
import pytest
from matplotlib import rc
from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.boltztrap import BoltztrapAnalyzer
from pymatgen.electronic_structure.cohp import CompleteCohp
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.plotter import (
    BoltztrapPlotter,
    BSDOSPlotter,
    BSPlotter,
    BSPlotterProjected,
    CohpPlotter,
    DosPlotter,
    fold_point,
    plot_brillouin_zone,
    plot_ellipsoid,
)
from pymatgen.io.vasp import Vasprun
from pymatgen.util.testing import TEST_FILES_DIR, VASP_IN_DIR, VASP_OUT_DIR, MatSciTest

BAND_TEST_DIR = f"{TEST_FILES_DIR}/electronic_structure/bandstructure"

rc("text", usetex=False)  # Disabling latex is needed for this test to work.


class TestDosPlotter(MatSciTest):
    def setup_method(self):
        with open(f"{BAND_TEST_DIR}/../dos/complete_dos.json", "rb") as file:
            self.dos = CompleteDos.from_dict(orjson.loads(file.read()))
            self.plotter = DosPlotter(sigma=0.2, stack=True)

    def test_add_dos_dict(self):
        dct = self.plotter.get_dos_dict()
        assert len(dct) == 0
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        dct = self.plotter.get_dos_dict()
        assert len(dct) == 4

    def test_get_dos_dict(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        dct = self.plotter.get_dos_dict()
        assert list(dct) == ["Li", "Fe", "P", "O"]

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_get_plot(self):
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        ax = self.plotter.get_plot()
        assert len(ax.lines) == 2
        out_path = f"{self.tmp_path}/dosplot.png"
        self.plotter.save_plot(out_path)
        assert os.path.isfile(out_path)

    def test_get_plot_limits(self):
        # Tests limit determination and if inverted_axes case
        # reproduces the same energy and DOS axis limits
        self.plotter.add_dos_dict(self.dos.get_element_dos(), key_sort_func=lambda x: x.X)
        # Contains energy and DOS limits and expected results
        with open(f"{BAND_TEST_DIR}/../plotter/complete_dos_limits.json", "rb") as file:
            limits_results = orjson.loads(file.read())

        for item in limits_results:
            ax = self.plotter.get_plot(xlim=item["energy_limit"], ylim=item["DOS_limit"])
            param_dict = self.get_plot_attributes(ax)
            ax_invert = self.plotter.get_plot(invert_axes=True, xlim=item["DOS_limit"], ylim=item["energy_limit"])
            param_dict_invert = self.get_plot_attributes(ax_invert)
            assert item["energy_result"] == approx(param_dict["xaxis_limits"])
            assert item["energy_result"] == approx(param_dict_invert["yaxis_limits"])
            assert item["DOS_result"] == approx(param_dict["yaxis_limits"])
            assert item["DOS_result"] == approx(param_dict_invert["xaxis_limits"])

    @staticmethod
    def get_plot_attributes(ax: plt.Axes):
        return {
            "xaxis_limits": list(ax.get_xlim()),
            "yaxis_limits": list(ax.get_ylim()),
        }


class TestBSPlotter(MatSciTest):
    def setup_method(self):
        with open(f"{BAND_TEST_DIR}/CaO_2605_bandstructure.json", encoding="utf-8") as file:
            dct = orjson.loads(file.read())
            self.bs = BandStructureSymmLine.from_dict(dct)
            self.plotter = BSPlotter(self.bs)

        assert len(self.plotter._bs) == 1, "wrong number of band objects"

        with open(f"{BAND_TEST_DIR}/N2_12103_bandstructure.json", encoding="utf-8") as file:
            dct = orjson.loads(file.read())
            self.sbs_sc = BandStructureSymmLine.from_dict(dct)

        with open(f"{BAND_TEST_DIR}/C_48_bandstructure.json", encoding="utf-8") as file:
            dct = orjson.loads(file.read())
            self.sbs_met = BandStructureSymmLine.from_dict(dct)

        self.plotter_multi = BSPlotter([self.sbs_sc, self.sbs_met])
        assert len(self.plotter_multi._bs) == 2, "wrong number of band objects"
        assert self.plotter_multi._nb_bands == [96, 96], "wrong number of bands"

    def test_add_bs(self):
        self.plotter_multi.add_bs(self.sbs_sc)
        assert len(self.plotter_multi._bs) == 3, "wrong number of band objects"
        assert self.plotter_multi._nb_bands == [96, 96, 96], "wrong number of bands"

    def test_get_branch_steps(self):
        steps_idx = BSPlotter._get_branch_steps(self.sbs_sc.branches)
        assert steps_idx == [0, 121, 132, 143], "wrong list of steps idx"

    def test_rescale_distances(self):
        rescaled_distances = self.plotter_multi._rescale_distances(self.sbs_sc, self.sbs_met)
        assert len(rescaled_distances) == len(self.sbs_met.distance), "wrong length of distances list"
        assert rescaled_distances[-1] == approx(6.5191398067252875), "wrong last distance value"
        assert rescaled_distances[148] == self.sbs_sc.distance[19], "wrong distance at high symm k-point"

    def test_interpolate_bands(self):
        data = self.plotter.bs_plot_data()
        dct = data["distances"]
        en = data["energy"]["1"]
        int_distances, int_energies = self.plotter._interpolate_bands(dct, en)

        assert len(int_distances) == 10, "wrong length of distances list"
        assert len(int_distances[0]) == 100, "wrong length of distances in a branch"
        assert len(int_energies) == 10, "wrong length of distances list"
        assert int_energies[0].shape == (16, 100), "wrong length of distances list"

    def test_bs_plot_data(self):
        assert len(self.plotter.bs_plot_data()["distances"]) == 10, "wrong number of sequences of branches"
        assert len(self.plotter.bs_plot_data()["distances"][0]) == 16, (
            "wrong number of distances in the first sequence of branches"
        )
        assert sum(len(dist) for dist in self.plotter.bs_plot_data()["distances"]) == 160, "wrong number of distances"

        length = len(self.plotter.bs_plot_data(split_branches=False)["distances"][0])
        assert length == 144, "wrong number of distances in the first sequence of branches"

        length = len(self.plotter.bs_plot_data(split_branches=False)["distances"])
        assert length == 2, "wrong number of distances in the first sequence of branches"

        assert self.plotter.bs_plot_data()["ticks"]["label"][5] == "K", "wrong tick label"
        assert len(self.plotter.bs_plot_data()["ticks"]["label"]) == 19, "wrong number of tick labels"

    def test_get_ticks(self):
        assert self.plotter.get_ticks()["label"][5] == "K", "wrong tick label"
        assert self.plotter.get_ticks()["distance"][5] == approx(2.406607625322699), "wrong tick distance"

    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_get_plot(self):
        # zero_to_efermi = True, ylim = None, smooth = False,
        # vbm_cbm_marker = False, smooth_tol = None

        ax = self.plotter.get_plot()
        assert_allclose(ax.get_ylim(), (-4.0, 7.6348)), "wrong ylim"
        ax = self.plotter.get_plot(smooth=True)
        ax = self.plotter.get_plot(vbm_cbm_marker=True)
        self.plotter.save_plot(f"{self.tmp_path}/bsplot.png")
        assert os.path.isfile(f"{self.tmp_path}/bsplot.png")
        plt.close("all")

        # test plotter with 2 bandstructures
        ax = self.plotter_multi.get_plot()
        assert len(ax.get_lines()) == 874, "wrong number of lines"
        assert_allclose(ax.get_ylim(), (-10.0, 10.0)), "wrong ylim"
        ax = self.plotter_multi.get_plot(zero_to_efermi=False)
        assert_allclose(ax.get_ylim(), (-15.2379, 12.67141266)), "wrong ylim"
        ax = self.plotter_multi.get_plot(smooth=True)
        self.plotter_multi.save_plot(f"{self.tmp_path}/bsplot.png")
        assert os.path.isfile(f"{self.tmp_path}/bsplot.png")
        plt.close("all")


class TestBSPlotterProjected:
    def setup_method(self):
        with open(f"{BAND_TEST_DIR}/Cu2O_361_bandstructure.json", "rb") as file:
            self.bs_Cu2O = BandStructureSymmLine.from_dict(orjson.loads(file.read()))
        self.plotter_Cu2O = BSPlotterProjected(self.bs_Cu2O)

        with open(f"{TEST_FILES_DIR}/electronic_structure/boltztrap2/PbTe_bandstructure.json", "rb") as file:
            self.bs_PbTe = BandStructureSymmLine.from_dict(orjson.loads(file.read()))

    def test_methods(self):
        # Minimal baseline testing for get_plot. not a true test. Just checks that
        # it can actually execute.
        self.plotter_Cu2O.get_elt_projected_plots()
        self.plotter_Cu2O.get_elt_projected_plots_color()
        axes = self.plotter_Cu2O.get_projected_plots_dots({"Cu": ["d", "s"], "O": ["p"]})
        assert len(axes) == 4, f"{len(axes)=}"
        assert len(axes[0].get_lines()) == 4903, f"{len(axes[0].get_lines())=}"
        assert len(axes[-1].get_lines()) == 0, f"{len(axes[-1].get_lines())=}"
        ax = self.plotter_Cu2O.get_projected_plots_dots_patom_pmorb(
            {"Cu": ["dxy", "s", "px"], "O": ["px", "py", "pz"]},
            {"Cu": [3, 5], "O": [1]},
        )
        assert isinstance(ax, plt.Axes)
        assert len(ax.get_lines()) == 44_127
        assert ax.get_ylim() == approx((-4.0, 4.5047))

        with pytest.raises(
            ValueError,
            match="Can't plot projections on a band structure without projections data",
        ):
            BSPlotterProjected(self.bs_PbTe)


class TestBSDOSPlotter:
    # Minimal baseline testing for get_plot. not a true test. Just checks that
    # it can actually execute.
    def test_methods(self):
        vasp_run = Vasprun(f"{VASP_OUT_DIR}/vasprun_Si_bands.xml.gz")
        plotter = BSDOSPlotter()
        band_struct = vasp_run.get_band_structure(kpoints_filename=f"{VASP_IN_DIR}/KPOINTS_Si_bands")
        ax = plotter.get_plot(band_struct)
        assert isinstance(ax, plt.Axes)
        plt.close()
        bs_ax, dos_ax = plotter.get_plot(band_struct, vasp_run.complete_dos)
        assert isinstance(bs_ax, plt.Axes)
        assert isinstance(dos_ax, plt.Axes)
        plt.close("all")

        with open(f"{TEST_FILES_DIR}/electronic_structure/plotter/SrBa2Sn2O7.json", "rb") as file:
            band_struct_dict = orjson.loads(file.read())
        # generate random projections
        data_structure = [[[[0 for _ in range(12)] for _ in range(9)] for _ in range(70)] for _ in range(90)]
        band_struct_dict["projections"]["1"] = data_structure
        dct = band_struct_dict["projections"]["1"]
        rng = np.random.default_rng()
        for ii in range(len(dct)):
            for jj in range(len(dct[ii])):
                for kk in range(len(dct[ii][jj])):
                    for ll in range(len(dct[ii][jj][kk])):
                        dct[ii][jj][kk][ll] = 0
                        # d[i][j][k][m] = rng.random()
                    # generate random number for two atoms
                    a = rng.integers(0, 7)
                    b = rng.integers(0, 7)
                    # c = rng.integers(0, 7)
                    dct[ii][jj][kk][a] = rng.random()
                    dct[ii][jj][kk][b] = rng.random()
                    # d[i][j][k][c] = rng.random()
        band_struct = BandStructureSymmLine.from_dict(band_struct_dict)
        ax = plotter.get_plot(band_struct)
        assert isinstance(ax, plt.Axes)


class TestPlotBZ:
    def setup_method(self):
        self.rec_latt = Structure.from_file(f"{TEST_FILES_DIR}/io/cssr/Si.cssr").lattice.reciprocal_lattice
        self.kpath = [[[0.0, 0.0, 0.0], [0.5, 0.0, 0.5], [0.5, 0.25, 0.75], [0.375, 0.375, 0.75]]]
        self.labels = {
            "\\Gamma": [0.0, 0.0, 0.0],
            "K": [0.375, 0.375, 0.75],
            "L": [0.5, 0.5, 0.5],
            "U": [0.625, 0.25, 0.625],
            "W": [0.5, 0.25, 0.75],
            "X": [0.5, 0.0, 0.5],
        }
        self.hessian = [
            [17.64757034, 3.90159625, -4.77845607],
            [3.90159625, 14.88874142, 6.75776076],
            [-4.77845607, 6.75776076, 12.12987493],
        ]
        self.center = [0.41, 0.0, 0.41]
        self.points = [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]

    def test_bz_plot(self):
        _, ax = plot_ellipsoid(self.hessian, self.center, lattice=self.rec_latt)
        plot_brillouin_zone(
            self.rec_latt,
            lines=self.kpath,
            labels=self.labels,
            kpoints=self.points,
            ax=ax,
            show=False,
        )

    def test_fold_point(self):
        assert_allclose(
            fold_point([0.0, -0.5, 0.5], lattice=self.rec_latt),
            self.rec_latt.get_cartesian_coords([0.0, 0.5, 0.5]),
            atol=1e-8,
        )
        assert_allclose(
            fold_point([0.1, -0.6, 0.2], lattice=self.rec_latt),
            self.rec_latt.get_cartesian_coords([0.1, 0.4, 0.2]),
            atol=1e-8,
        )


@pytest.mark.xfail(reason="TODO: need someone to fix this")
@pytest.mark.skipif(not which("x_trans"), reason="No x_trans executable found")
class TestBoltztrapPlotter:
    def setup_method(self):
        bz = BoltztrapAnalyzer.from_files(f"{TEST_FILES_DIR}/boltztrap/transp/")
        self.plotter = BoltztrapPlotter(bz)

    def test_plot_carriers(self):
        ax = self.plotter.plot_carriers()
        assert len(ax.get_lines()) == 7, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.0702422655947665), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(6.525490122298364e22), "wrong 1 data in line 0"
        plt.close()

    def test_plot_complexity_factor_mu(self):
        pytest.importorskip("fdint")
        ax = self.plotter.plot_complexity_factor_mu()
        assert len(ax.get_lines()) == 2, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.0702422655947665), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(0.004708835456903449), "wrong 1 data in line 0"
        plt.close()

    def test_plot_conductivity_dop(self):
        ax = self.plotter.plot_conductivity_dop()
        assert len(ax.get_lines()) == 8, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(1e15), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(0.3801957596666667), "wrong 1 data in line 0"
        plt.close()

    def test_plot_conductivity_mu(self):
        ax = self.plotter.plot_conductivity_mu()
        assert len(ax.get_lines()) == 9, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.0702422655947665), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(1965.1306), "wrong 1 data in line 0"
        plt.close()

    def test_plot_conductivity_temp(self):
        ax = self.plotter.plot_conductivity_temp()
        assert len(ax.get_lines()) == 6, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(100), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(0.3801957596666667), "wrong 1 data in line 0"
        plt.close()

    def test_plot_dos(self):
        ax = self.plotter.plot_dos()
        assert len(ax.get_lines()) == 3, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.4197044934588674), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(0.0), "wrong 1 data in line 0"
        plt.close()

    def test_plot_eff_mass_dop(self):
        ax = self.plotter.plot_eff_mass_dop()
        assert len(ax.get_lines()) == 8, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(1e15), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(1.4231240011719886), "wrong 1 data in line 0"
        plt.close()

    def test_plot_eff_mass_temp(self):
        ax = self.plotter.plot_eff_mass_temp()
        assert len(ax.get_lines()) == 6, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(100), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(1.4231240011719886), "wrong 1 data in line 0"
        plt.close()

    def test_plot_hall_carriers(self):
        ax = self.plotter.plot_hall_carriers()
        assert len(ax.get_lines()) == 7, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.0702422655947665), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(9.538187273102463e17), "wrong 1 data in line 0"
        plt.close()

    def test_plot_power_factor_dop(self):
        ax = self.plotter.plot_power_factor_dop()
        assert len(ax.get_lines()) == 8, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(1e15), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(0.40606868935796925), "wrong 1 data in line 0"
        plt.close()

    def test_plot_power_factor_mu(self):
        ax = self.plotter.plot_power_factor_mu()
        assert len(ax.get_lines()) == 9, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.0702422655947665), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(365.5514594136157), "wrong 1 data in line 0"
        plt.close()

    def test_plot_power_factor_temp(self):
        ax = self.plotter.plot_power_factor_temp()
        assert len(ax.get_lines()) == 6, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(100), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(0.40606868935796925), "wrong 1 data in line 0"
        plt.close()

    def test_plot_seebeck_dop(self):
        ax = self.plotter.plot_seebeck_dop()
        assert len(ax.get_lines()) == 8, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(1e15), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(1050.8197666666667), "wrong 1 data in line 0"
        plt.close()

    def test_plot_seebeck_eff_mass_mu(self):
        pytest.importorskip("fdint")
        ax = self.plotter.plot_seebeck_eff_mass_mu()
        assert len(ax.get_lines()) == 2, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.0702422655947665), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(6412.881888198197), "wrong 1 data in line 0"
        plt.close()

    def test_plot_seebeck_mu(self):
        ax = self.plotter.plot_seebeck_mu()
        assert len(ax.get_lines()) == 9, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.0702422655947665), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(-433.11096000000003), "wrong 1 data in line 0"
        plt.close()

    def test_plot_seebeck_temp(self):
        ax = self.plotter.plot_seebeck_temp()
        assert len(ax.get_lines()) == 6, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(100), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(1050.8197666666667), "wrong 1 data in line 0"
        plt.close()

    def test_plot_zt_dop(self):
        ax = self.plotter.plot_zt_dop()
        assert len(ax.get_lines()) == 8, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(1e15), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(4.060682863129955e-05), "wrong 1 data in line 0"
        plt.close()

    def test_plot_zt_mu(self):
        ax = self.plotter.plot_zt_mu()
        assert len(ax.get_lines()) == 9, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(-2.0702422655947665), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(0.2153839699235254), "wrong 1 data in line 0"
        plt.close()

    def test_plot_zt_temp(self):
        ax = self.plotter.plot_zt_temp()
        assert len(ax.get_lines()) == 6, "wrong number of lines"
        assert ax.get_lines()[0].get_data()[0][0] == approx(100), "wrong 0 data in line 0"
        assert ax.get_lines()[0].get_data()[1][0] == approx(4.060682863129955e-05), "wrong 1 data in line 0"
        plt.close()


class TestCohpPlotter(MatSciTest):
    def setup_method(self):
        path = f"{TEST_FILES_DIR}/electronic_structure/cohp/complete_cohp_lobster.json"
        with open(path, "rb") as file:
            self.cohp = CompleteCohp.from_dict(orjson.loads(file.read()))
        path = f"{TEST_FILES_DIR}/electronic_structure/cohp/complete_coop_lobster.json"
        with open(path, "rb") as file:
            self.coop = CompleteCohp.from_dict(orjson.loads(file.read()))
        self.cohp_plot = CohpPlotter(zero_at_efermi=False)
        self.coop_plot = CohpPlotter(are_coops=True)

    def test_attributes(self):
        assert not self.cohp_plot.are_coops
        assert self.coop_plot.are_coops
        assert not self.cohp_plot.zero_at_efermi
        assert self.coop_plot.zero_at_efermi
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        cohp_energies = self.cohp_plot._cohps["1"]["energies"]
        assert len(cohp_energies) == 301
        assert cohp_energies[0] == approx(-0.27768)
        assert cohp_energies[-1] == approx(14.77248)
        self.coop_plot.add_cohp_dict(self.coop.all_cohps)
        coop_energies = self.coop_plot._cohps["10"]["energies"]
        assert len(coop_energies) == 241
        assert coop_energies[0] == approx(-6.02510)
        assert coop_energies[-1] == approx(6.02510)

    def test_add_cohp_dict(self):
        # Sorts the populations by z-coordinates of the sites
        def sortkeys(sites):
            return sites[0].z, sites[1].z

        sorted_keys = ["3", "4", "7", "8", "9", "10", "11", "6", "5", "2", "1"]

        d_coop = self.coop_plot.get_cohp_dict()
        assert len(d_coop) == 0
        bonds = self.coop.bonds
        self.coop_plot.add_cohp_dict(self.coop.all_cohps, key_sort_func=lambda x: sortkeys(bonds[x]["sites"]))
        d_coop = self.coop_plot.get_cohp_dict()
        assert len(d_coop) == 11
        assert list(self.coop_plot._cohps) == sorted_keys

    def test_get_cohp_dict(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        d_cohp = self.cohp_plot.get_cohp_dict()
        for bond in ["1", "2"]:
            assert bond in d_cohp

    def test_get_plot(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        ax_cohp = self.cohp_plot.get_plot()

        assert ax_cohp.get_xlabel() == "-COHP"
        assert ax_cohp.get_ylabel() == "$E$ (eV)"
        legend_labels = ax_cohp.get_legend_handles_labels()[1]
        assert len(self.cohp_plot._cohps) == len(legend_labels)
        assert ax_cohp.lines[0].get_linestyle() == "-"
        assert ax_cohp.lines[1].get_linestyle() == "--"
        for label in legend_labels:
            assert label in self.cohp_plot._cohps
        lines_index = legend_labels.index("1")
        line_styles = {Spin.up: "-", Spin.down: "--"}
        cohp_fe_fe = self.cohp.all_cohps["1"]
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * lines_index + s]
            assert_allclose(lines.get_xdata(), -cohp_fe_fe.cohp[spin])
            assert_allclose(lines.get_ydata(), self.cohp.energies)
            assert lines.get_linestyle() == line_styles[spin]
        plt.close()

        ax_cohp = self.cohp_plot.get_plot(invert_axes=False, plot_negative=False)

        assert ax_cohp.get_xlabel() == "$E$ (eV)"
        assert ax_cohp.get_ylabel() == "COHP"
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * lines_index + s]
            assert_allclose(lines.get_xdata(), self.cohp.energies)
            assert_allclose(lines.get_ydata(), cohp_fe_fe.cohp[spin])
        plt.close()

        ax_cohp = self.cohp_plot.get_plot(integrated=True)

        assert ax_cohp.get_xlabel() == "-ICOHP (eV)"
        for s, spin in enumerate([Spin.up, Spin.down]):
            lines = ax_cohp.lines[2 * lines_index + s]
            assert_allclose(lines.get_xdata(), -cohp_fe_fe.icohp[spin])

        coop_dict = {"Bi5-Bi6": self.coop.all_cohps["10"]}
        self.coop_plot.add_cohp_dict(coop_dict)
        ax_coop = self.coop_plot.get_plot()
        assert ax_coop.get_xlabel() == "COOP"
        assert ax_coop.get_ylabel() == "$E - E_f$ (eV)"
        lines_coop = ax_coop.get_lines()[0]
        assert_allclose(lines_coop.get_ydata(), self.coop.energies - self.coop.efermi)
        coop_bi_bi = self.coop.all_cohps["10"].cohp[Spin.up]
        assert_allclose(lines_coop.get_xdata(), coop_bi_bi)

        # cleanup
        plt.close("all")

    def test_save_plot(self):
        self.cohp_plot.add_cohp_dict(self.cohp.all_cohps)
        ax = self.cohp_plot.get_plot()
        assert isinstance(ax, plt.Axes)
        self.cohp_plot.save_plot(f"{self.tmp_path}/cohpplot.png")
        assert os.path.isfile(f"{self.tmp_path}/cohpplot.png")
        plt.close("all")
