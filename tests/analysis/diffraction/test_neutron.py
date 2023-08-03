from __future__ import annotations

import pytest
from pytest import approx

from pymatgen.analysis.diffraction.neutron import NDCalculator
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

"""
These calculated values were verified with VESTA and FullProf.
"""

__author__ = "Yuta Suzuki"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Yuta Suzuki"
__email__ = "resnant@outlook.jp"
__date__ = "4/19/18"


class TestNDCalculator(PymatgenTest):
    def test_get_pattern(self):
        struct = self.get_structure("CsCl")
        c = NDCalculator(wavelength=1.54184)  # CuKa radiation
        nd = c.get_pattern(struct, two_theta_range=(0, 90))
        # Check the first two peaks
        assert nd.x[0] == approx(21.107738329639844)
        assert nd.hkls[0] == [{"hkl": (1, 0, 0), "multiplicity": 6}]
        assert nd.d_hkls[0] == approx(4.2089999999999996)
        assert nd.x[1] == approx(30.024695921112777)
        assert nd.hkls[1] == [{"hkl": (1, 1, 0), "multiplicity": 12}]
        assert nd.d_hkls[1] == approx(2.976212442014178)

        struct = self.get_structure("LiFePO4")
        nd = c.get_pattern(struct, two_theta_range=(0, 90))
        assert nd.x[1] == approx(17.03504233621785)
        assert nd.y[1] == approx(46.2985965)

        struct = self.get_structure("Li10GeP2S12")
        nd = c.get_pattern(struct, two_theta_range=(0, 90))
        assert nd.x[1] == approx(14.058274883353876)
        assert nd.y[1] == approx(3.60588013)

        # Test a hexagonal structure.
        struct = self.get_structure("Graphite")
        nd = c.get_pattern(struct, two_theta_range=(0, 90))
        assert nd.x[0] == approx(26.21057350859598)
        assert nd.y[0] == approx(100)
        assert nd.x[2] == approx(44.39599754)
        assert nd.y[2] == approx(42.62382267)
        assert len(nd.hkls[0][0]) == approx(2)

        # Test an exception in case of the input element is
        # not in scattering length table.
        # This curium structure is just for test, not the actual structure.
        something = Structure(Lattice.cubic(a=1), ["Cm"], [[0, 0, 0]])
        with pytest.raises(
            ValueError, match="Unable to calculate ND pattern as there is no scattering coefficients for Cm."
        ):
            nd = c.get_pattern(something, two_theta_range=(0, 90))

        # Test with Debye-Waller factor
        struct = self.get_structure("Graphite")
        c = NDCalculator(wavelength=1.54184, debye_waller_factors={"C": 1})
        nd = c.get_pattern(struct, two_theta_range=(0, 90))
        assert nd.x[0] == approx(26.21057350859598)
        assert nd.y[0] == approx(100)
        assert nd.x[2] == approx(44.39599754)
        assert nd.y[2] == approx(39.471514740)

    def test_get_plot(self):
        struct = self.get_structure("Graphite")
        c = NDCalculator(wavelength=1.54184, debye_waller_factors={"C": 1})
        c.get_plot(struct, two_theta_range=(0, 90))
