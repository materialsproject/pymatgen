from __future__ import annotations

import pytest
from pytest import approx

from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.util.testing import PymatgenTest

"""
TODO: Modify unittest doc.
"""

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "5/22/14"


class TestXRDCalculator(PymatgenTest):
    def test_type_wavelength(self):
        """Test TypeError is raised if wavelength is unaccepted type."""
        wavelength = [1.78, 2.78]  # just a list
        with pytest.raises(TypeError) as exc:
            XRDCalculator(wavelength)

        assert "type(wavelength)=<class 'list'> must be either float, int or str" in str(exc.value)

    def test_get_pattern(self):
        struct = self.get_structure("CsCl")
        xrd_calc = XRDCalculator()
        xrd = xrd_calc.get_pattern(struct, two_theta_range=(0, 90))
        assert xrd.to_json()  # Test MSONable property
        # Check the first two peaks
        assert xrd.x[0] == approx(21.107738329639844)
        assert xrd.y[0] == approx(36.483184003748946)
        assert xrd.hkls[0] == [{"hkl": (1, 0, 0), "multiplicity": 6}]
        assert xrd.d_hkls[0] == approx(4.2089999999999996)
        assert xrd.x[1] == approx(30.024695921112777)
        assert xrd.y[1] == approx(100)
        assert xrd.hkls[1] == [{"hkl": (1, 1, 0), "multiplicity": 12}]
        assert xrd.d_hkls[1] == approx(2.976212442014178)

        struct = self.get_structure("LiFePO4")
        xrd = xrd_calc.get_pattern(struct, two_theta_range=(0, 90))
        assert xrd.x[1] == approx(17.03504233621785)
        assert xrd.y[1] == approx(50.400928948337075)

        struct = self.get_structure("Li10GeP2S12")
        xrd = xrd_calc.get_pattern(struct, two_theta_range=(0, 90))
        assert xrd.x[1] == approx(14.058274883353876)
        assert xrd.y[1] == approx(4.4111123641667671)

        # Test a hexagonal structure.
        struct = self.get_structure("Graphite")

        xrd = xrd_calc.get_pattern(struct, two_theta_range=(0, 90))
        assert xrd.x[0] == approx(26.21057350859598)
        assert xrd.y[0] == approx(100)
        assert len(xrd.hkls[0][0]["hkl"]) == 4

        # Add test case with different lengths of coefficients.
        # Also test d_hkl.
        coords = [
            [0.25, 0.25, 0.173],
            [0.75, 0.75, 0.827],
            [0.75, 0.25, 0],
            [0.25, 0.75, 0],
            [0.25, 0.25, 0.676],
            [0.75, 0.75, 0.324],
        ]
        sp = ["Si", "Si", "Ru", "Ru", "Pr", "Pr"]
        struct = Structure(Lattice.tetragonal(4.192, 6.88), sp, coords)
        xrd = xrd_calc.get_pattern(struct)
        assert xrd.x[0] == approx(12.86727341476735)
        assert xrd.y[0] == approx(31.448239816769796)
        assert xrd.d_hkls[0] == approx(6.88)
        assert len(xrd) == 42
        xrd = xrd_calc.get_pattern(struct, two_theta_range=[0, 60])
        assert len(xrd) == 18

        # Test with and without Debye-Waller factor
        tungsten = Structure(Lattice.cubic(3.1653), ["W"] * 2, [[0, 0, 0], [0.5, 0.5, 0.5]])
        xrd = xrd_calc.get_pattern(tungsten, scaled=False)
        assert xrd.x[0] == approx(40.294828554672264)
        assert xrd.y[0] == approx(2414237.5633093244)
        assert xrd.d_hkls[0] == approx(2.2382050944897789)
        xrd_calc = XRDCalculator(debye_waller_factors={"W": 0.1526})
        xrd = xrd_calc.get_pattern(tungsten, scaled=False)
        assert xrd.x[0] == approx(40.294828554672264)
        assert xrd.y[0] == approx(2377745.2296686019)
        assert xrd.d_hkls[0] == approx(2.2382050944897789)
