from __future__ import annotations

import math
import shutil

import numpy as np
import orjson
import pytest
from matplotlib import pyplot as plt
from numpy.testing import assert_allclose

from pymatgen.analysis.transition_state import NEBAnalysis, combine_neb_plots
from pymatgen.io.vasp import Outcar
from pymatgen.util.testing import MatSciTest
from tests.testing import TEST_FILES_DIR

"""
TODO: Modify unittest doc.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__date__ = "2/5/16"


TEST_DIR = f"{TEST_FILES_DIR}/io/vasp/fixtures/neb_analysis"


class TestNEBAnalysis(MatSciTest):
    def test_run(self):
        neb_analysis1 = NEBAnalysis.from_dir(f"{TEST_DIR}/neb1/neb")
        neb_analysis1_from_dict = NEBAnalysis.from_dict(neb_analysis1.as_dict())
        json_data = orjson.dumps(neb_analysis1.as_dict()).decode()

        neb_dict = orjson.loads(json_data)
        neb_analysis1_from_json_data = NEBAnalysis.from_dict(neb_dict)

        assert_allclose(neb_analysis1.energies[0], -255.97992669000001)
        assert_allclose(neb_analysis1.energies[3], -255.84261996000001)
        assert_allclose(neb_analysis1.r, neb_analysis1_from_dict.r)
        assert_allclose(neb_analysis1.energies, neb_analysis1_from_dict.energies)
        assert_allclose(neb_analysis1.forces, neb_analysis1_from_dict.forces)
        assert neb_analysis1.structures == neb_analysis1_from_dict.structures

        assert_allclose(neb_analysis1.r, neb_analysis1_from_json_data.r)
        assert_allclose(neb_analysis1.energies, neb_analysis1_from_json_data.energies)
        assert_allclose(neb_analysis1.forces, neb_analysis1_from_json_data.forces)
        assert neb_analysis1.structures == neb_analysis1_from_json_data.structures

        assert_allclose(neb_analysis1.get_extrema()[1][0], (0.50023335723480078, 325.20043063935128))

        neb_analysis1.setup_spline(zero_slope_saddle=True)
        assert_allclose(neb_analysis1.get_extrema()[1][0], (0.50023335723480078, 325.20003984140203))
        with open(f"{TEST_DIR}/neb2/neb_analysis2.json", "rb") as file:
            neb_analysis2_dict = orjson.loads(file.read())
        neb_analysis2 = NEBAnalysis.from_dict(neb_analysis2_dict)
        assert_allclose(neb_analysis2.get_extrema()[1][0], (0.37255257367467326, 562.40825334519991))

        neb_analysis2.setup_spline(zero_slope_saddle=True)
        assert_allclose(neb_analysis2.get_extrema()[1][0], (0.30371133723478794, 528.46229631648691))

    def test_combine_neb_plots(self):
        neb_dir = f"{TEST_DIR}/neb1/neb"
        neb_analysis = NEBAnalysis.from_dir(neb_dir)
        combine_neb_plots([neb_analysis, neb_analysis])

    def test_get_plot(self):
        neb_dir = f"{TEST_DIR}/neb1/neb"
        neb_analysis = NEBAnalysis.from_dir(neb_dir)
        ax = neb_analysis.get_plot()
        assert isinstance(ax, plt.Axes)
        assert math.isclose(neb_analysis.get_extrema()[1][0][1], 325.20043063935128)
        assert ax.texts[0].get_text() == "325 meV", "Unexpected annotation text"
        assert ax.get_xlabel() == "Reaction Coordinate"
        assert ax.get_ylabel() == "Energy (meV)"


NEB1_DIR = TEST_FILES_DIR / "io" / "vasp" / "fixtures" / "neb_analysis" / "neb1" / "neb"
# Energies (eV) of images 00-04 of the neb1 run
NEB1_ENERGIES = [-255.97992669, -255.84212624, -255.65472488, -255.84261996, -255.98097765]


@pytest.fixture(scope="module")
def _neb1_parsed() -> NEBAnalysis:
    return NEBAnalysis.from_dir(NEB1_DIR)


@pytest.fixture
def neb1(_neb1_parsed) -> NEBAnalysis:
    """Fresh copy of the parsed neb1 run, so tests may modify it."""
    neb = _neb1_parsed
    return NEBAnalysis(neb.r.copy(), neb.energies.copy(), neb.forces.copy(), list(neb.structures))


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _shifted(neb: NEBAnalysis, shift: float = 0, reverse: bool = False) -> NEBAnalysis:
    """Copy of an NEBAnalysis with energies shifted by `shift` eV, optionally traversed backwards."""
    energies = np.asarray(neb.energies) + shift
    if not reverse:
        return NEBAnalysis(neb.r, energies, neb.forces, list(neb.structures))
    return NEBAnalysis(neb.r[-1] - neb.r[::-1], energies[::-1], neb.forces[::-1], neb.structures[::-1])


class TestNEBAnalysisExtra:
    def test_from_dir_values(self, neb1):
        assert_allclose(neb1.energies, NEB1_ENERGIES)
        assert_allclose(neb1.relative_energies, np.array(NEB1_ENERGIES) - NEB1_ENERGIES[0])
        # Terminal images carry no tangent force
        assert_allclose(neb1.forces, [0, -0.539991, 0.007216, 0.540828, 0])
        assert_allclose(neb1.r, [0, 0.74841907, 1.49822886, 2.24944916, 2.99860051], rtol=1e-7)
        assert len(neb1.structures) == 5
        assert neb1.zero_slope_saddle is False

    @pytest.mark.parametrize("subdirs", [("start", "end"), ("initial", "final")])
    def test_from_dir_terminal_subdirs(self, tmp_path, subdirs):
        root = tmp_path / "neb"
        shutil.copytree(NEB1_DIR, root)
        for src, sub in zip(("00", "04"), subdirs, strict=True):
            (root / sub).mkdir()
            shutil.move(root / src / "OUTCAR", root / sub / "OUTCAR")
        neb = NEBAnalysis.from_dir(root)
        assert_allclose(neb.energies, NEB1_ENERGIES)
        assert len(neb.structures) == 5

    def test_from_dir_relaxation_dirs(self, tmp_path, neb1):
        root = tmp_path / "neb"
        shutil.copytree(NEB1_DIR, root)
        # Explicit relaxation_dirs take precedence over the OUTCARs in 00 and 04
        relaxation_dirs = (tmp_path / "relax_start", tmp_path / "relax_end")
        for src, target in zip(("01", "03"), relaxation_dirs, strict=True):
            target.mkdir()
            shutil.copy(root / src / "OUTCAR", target / "OUTCAR")
        neb = NEBAnalysis.from_dir(root, relaxation_dirs=relaxation_dirs)
        assert_allclose(neb.energies, [NEB1_ENERGIES[1], *NEB1_ENERGIES[1:4], NEB1_ENERGIES[3]])
        # Terminal structures are still read from the POSCARs in 00 and 04
        assert neb.structures == neb1.structures

    def test_from_dir_missing_terminal_outcar(self, tmp_path):
        root = tmp_path / "neb"
        shutil.copytree(NEB1_DIR, root)
        (root / "04" / "OUTCAR").unlink()
        with pytest.raises(ValueError, match="OUTCAR cannot be found for terminal point"):
            NEBAnalysis.from_dir(root)

    def test_from_outcars(self, neb1):
        outcars = [Outcar(NEB1_DIR / f"0{idx}" / "OUTCAR") for idx in range(5)]
        neb = NEBAnalysis.from_outcars(outcars, neb1.structures, zero_slope_saddle=True)
        assert neb.zero_slope_saddle
        assert_allclose(neb.energies, NEB1_ENERGIES)
        assert_allclose(neb.r, neb1.r)
        with pytest.raises(ValueError, match="# of Outcars must be same as # of Structures"):
            NEBAnalysis.from_outcars(outcars, neb1.structures[:-1])

    def test_spline_options_deprecated(self, neb1):
        with pytest.warns(DeprecationWarning, match="spline_options"):
            neb = NEBAnalysis(
                neb1.r, neb1.energies, neb1.forces, neb1.structures, spline_options={"saddle_point": "zero_slope"}
            )
        assert neb.zero_slope_saddle
        assert_allclose(neb.get_extrema()[1][0], (0.50023335723480078, 325.20003984140203))

        with pytest.warns(DeprecationWarning, match="spline_options"):
            neb.setup_spline(spline_options={})
        assert neb.zero_slope_saddle is False
        assert_allclose(neb.get_extrema()[1][0], (0.50023335723480078, 325.20043063935128))

    @pytest.mark.parametrize("zero_slope_saddle", [False, True])
    def test_spline_passes_through_images(self, neb1, zero_slope_saddle):
        neb1.setup_spline(zero_slope_saddle=zero_slope_saddle)
        assert_allclose(neb1.spline(neb1.r), neb1.relative_energies, atol=1e-12)
        # Clamped (zero slope) end points
        assert_allclose(neb1.spline(neb1.r[[0, -1]], 1), 0, atol=1e-12)
        if zero_slope_saddle:
            assert neb1.spline(neb1.r[2], 1) == pytest.approx(0, abs=1e-12)

    def test_get_extrema_with_local_minimum(self, neb1):
        # Double-humped path with an intermediate at ~100 meV
        neb = NEBAnalysis(np.arange(5.0), [0, 0.3, 0.1, 0.4, 0], np.zeros(5), neb1.structures)
        minima, maxima = neb.get_extrema()
        assert len(minima) == 1
        assert len(maxima) == 2
        assert_allclose(minima[0], (0.4875, 97.709375))
        assert_allclose(maxima, [(0.26, 300.9088), (0.7375, 402.37991071)])
        minima_abs, maxima_abs = neb.get_extrema(normalize_rxn_coordinate=False)
        assert_allclose(minima_abs[0], (1.95, 97.709375))
        assert_allclose([m[0] for m in maxima_abs], [1.04, 2.95])

    def test_get_plot_options(self, neb1):
        ax = neb1.get_plot(normalize_rxn_coordinate=False, label_barrier=False)
        assert len(ax.texts) == 0
        points, curve = ax.lines
        assert_allclose(points.get_xdata(), neb1.r)
        assert_allclose(points.get_ydata(), neb1.relative_energies * 1000)
        assert curve.get_xdata()[-1] == pytest.approx(2.99, abs=0.01)

        ax = neb1.get_plot()
        barrier_line = ax.lines[-1]
        assert barrier_line.get_linestyle() == "--"
        assert_allclose(barrier_line.get_ydata(), [325.20043, 325.20043], rtol=1e-6)
        assert_allclose(barrier_line.get_xdata()[1], 0.5002333572348008, rtol=1e-6)
        assert ax.texts[0].get_text() == "325 meV"

    def test_as_from_dict(self, neb1):
        dct = neb1.as_dict()
        assert dct["@class"] == "NEBAnalysis"
        assert_allclose(dct["energies"], NEB1_ENERGIES)
        new = NEBAnalysis.from_dict(dct)
        assert_allclose(new.get_extrema()[1], neb1.get_extrema()[1])

    @pytest.mark.xfail(
        reason="NEBAnalysis.as_dict does not serialize zero_slope_saddle, so the spline mode is lost",
        strict=True,
    )
    def test_as_from_dict_zero_slope_saddle(self, neb1):
        neb = NEBAnalysis(neb1.r, neb1.energies, neb1.forces, neb1.structures, zero_slope_saddle=True)
        assert NEBAnalysis.from_dict(neb.as_dict()).zero_slope_saddle


class TestCombineNEBPlots:
    @pytest.mark.parametrize("reverse", [False, True])
    def test_combine_at_start(self, neb1, reverse):
        # Paths sharing neb1's start (as start or end of neb2): combined path runs
        # neb1 backwards (or reversed neb2), then forwards, i.e. 04 -> 00 -> 04
        combined = combine_neb_plots([neb1, _shifted(neb1, reverse=reverse)])
        assert_allclose(combined.energies, [NEB1_ENERGIES[i] for i in (4, 3, 2, 1, 0, 1, 2, 3, 4)])
        assert len(combined.structures) == 9
        assert np.all(np.diff(combined.r) > 0)
        assert_allclose(combined.r[-1], 2 * neb1.r[-1])
        assert len(neb1.structures) == 5

    def test_combine_end_start(self, neb1):
        # neb2 starts where neb1 ends: neb1 forwards then neb2 forwards
        delta = NEB1_ENERGIES[-1] - NEB1_ENERGIES[0]
        combined = combine_neb_plots([neb1, _shifted(neb1, delta)])
        assert_allclose(combined.energies, [*NEB1_ENERGIES, *(np.array(NEB1_ENERGIES[1:]) + delta)])
        assert_allclose(combined.r, [*neb1.r, *(neb1.r[1:] + neb1.r[-1])])
        assert_allclose(combined.forces, [*neb1.forces, *neb1.forces[1:]])
        assert len(combined.structures) == 9
        # The input NEBAnalysis must not be mutated
        assert len(neb1.structures) == 5

    def test_combine_end_end(self, neb1):
        # neb2 ends where neb1 ends: neb1 forwards then neb2 backwards
        delta = NEB1_ENERGIES[-1] - NEB1_ENERGIES[0]
        neb2 = _shifted(neb1, delta, reverse=True)
        assert neb2.energies[-1] == pytest.approx(NEB1_ENERGIES[-1])
        combined = combine_neb_plots([neb1, neb2])
        assert_allclose(combined.energies, [*NEB1_ENERGIES, *(np.array(NEB1_ENERGIES[1:]) + delta)])
        assert np.all(np.diff(combined.r) > 0)
        assert_allclose(combined.r[-1], 2 * neb1.r[-1])
        assert len(combined.structures) == 9
        assert len(neb1.structures) == 5

    def test_combine_arranged_and_reversed(self, neb1):
        combined = combine_neb_plots([neb1, neb1], arranged_neb_analyses=True)
        # Joint energy is the mean of neb1's end and neb2's start
        joint = (NEB1_ENERGIES[-1] + NEB1_ENERGIES[0]) / 2
        assert_allclose(combined.energies, [*NEB1_ENERGIES[:-1], joint, *NEB1_ENERGIES[1:]])
        assert_allclose(combined.r, [*neb1.r, *(neb1.r[1:] + neb1.r[-1])])
        assert len(combined.structures) == 9
        assert len(neb1.structures) == 5

        reversed_combined = combine_neb_plots([neb1, neb1], arranged_neb_analyses=True, reverse_plot=True)
        assert_allclose(reversed_combined.energies, combined.energies[::-1])
        assert_allclose(reversed_combined.r, combined.r[-1] - combined.r[::-1], atol=1e-12)
        assert reversed_combined.structures == combined.structures[::-1]

    def test_combine_no_match(self, neb1):
        other = [struct.copy() for struct in neb1.structures]
        for struct in other:
            struct.replace_species(dict.fromkeys(struct.composition.elements, "Xe"))
        neb2 = NEBAnalysis(neb1.r, neb1.energies, neb1.forces, other)
        with pytest.raises(ValueError, match="no matched structures for connection"):
            combine_neb_plots([neb1, neb2])
