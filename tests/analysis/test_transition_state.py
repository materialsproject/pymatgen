from __future__ import annotations

import math

import orjson
from matplotlib import pyplot as plt
from numpy.testing import assert_allclose

from pymatgen.analysis.transition_state import NEBAnalysis, combine_neb_plots
from pymatgen.util.testing import TEST_FILES_DIR, MatSciTest

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
