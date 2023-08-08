from __future__ import annotations

import json
import warnings

import numpy as np

from pymatgen.analysis.transition_state import NEBAnalysis, combine_neb_plots
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest

"""
TODO: Modify unittest doc.
"""

__author__ = "David Waroquiers"
__copyright__ = "Copyright 2016, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyam Dwaraknath"
__email__ = "shyamd@lbl.gov"
__date__ = "2/5/16"


test_dir = f"{TEST_FILES_DIR}/neb_analysis"


class TestNEBAnalysis(PymatgenTest):
    def setUp(self):
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.simplefilter("default")

    def runTest(self):
        neb_analysis1 = NEBAnalysis.from_dir(f"{test_dir}/neb1/neb")
        neb_analysis1_from_dict = NEBAnalysis.from_dict(neb_analysis1.as_dict())
        json_data = json.dumps(neb_analysis1.as_dict())

        neb_dict = json.loads(json_data)
        neb_analysis1_from_json_data = NEBAnalysis.from_dict(neb_dict)

        assert np.allclose(neb_analysis1.energies[0], -255.97992669000001)
        assert np.allclose(neb_analysis1.energies[3], -255.84261996000001)
        assert np.allclose(neb_analysis1.r, neb_analysis1_from_dict.r)
        assert np.allclose(neb_analysis1.energies, neb_analysis1_from_dict.energies)
        assert np.allclose(neb_analysis1.forces, neb_analysis1_from_dict.forces)
        assert neb_analysis1.structures == neb_analysis1_from_dict.structures

        assert np.allclose(neb_analysis1.r, neb_analysis1_from_json_data.r)
        assert np.allclose(neb_analysis1.energies, neb_analysis1_from_json_data.energies)
        assert np.allclose(neb_analysis1.forces, neb_analysis1_from_json_data.forces)
        assert neb_analysis1.structures == neb_analysis1_from_json_data.structures

        assert np.allclose(neb_analysis1.get_extrema()[1][0], (0.50023335723480078, 325.20043063935128))

        neb_analysis1.setup_spline(spline_options={"saddle_point": "zero_slope"})
        assert np.allclose(neb_analysis1.get_extrema()[1][0], (0.50023335723480078, 325.20003984140203))
        with open(f"{test_dir}/neb2/neb_analysis2.json") as f:
            neb_analysis2_dict = json.load(f)
        neb_analysis2 = NEBAnalysis.from_dict(neb_analysis2_dict)
        assert np.allclose(neb_analysis2.get_extrema()[1][0], (0.37255257367467326, 562.40825334519991))

        neb_analysis2.setup_spline(spline_options={"saddle_point": "zero_slope"})
        assert np.allclose(neb_analysis2.get_extrema()[1][0], (0.30371133723478794, 528.46229631648691))

    def test_combine_neb_plots(self):
        neb_dir = f"{test_dir}/neb1/neb"
        neb_analysis = NEBAnalysis.from_dir(neb_dir)
        combine_neb_plots([neb_analysis, neb_analysis])
