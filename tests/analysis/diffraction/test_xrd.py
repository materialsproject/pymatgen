from __future__ import annotations

import pytest
from pytest import approx

from pymatgen.analysis.diffraction.xrd import XRDCalculator
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Species
from pymatgen.core.structure import Structure
from pymatgen.util.testing import MatSciTest

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "5/22/14"


class TestXRDCalculator(MatSciTest):
    def test_type_wavelength(self):
        """Test TypeError is raised if wavelength is unaccepted type."""
        wavelength = [1.78, 2.78]  # just a list
        with pytest.raises(TypeError, match="must be either float, int or str"):
            XRDCalculator(wavelength)

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

    def test_get_pattern_disordered_species_dw(self):
        """Disordered fcc Cu-Au with oxidation-state species, partial occupancies
        and Debye-Waller factors: exercises the occupancy-weighted sum over
        site.species and the symbol-keyed scattering/DW lookups. Reference
        values were generated with the per-atom (pre grouped-by-element)
        implementation and must not change.
        """
        struct = Structure(
            Lattice.cubic(3.677),
            [{Species("Cu2+"): 0.5, Species("Au3+"): 0.5}] * 4,
            [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        )
        xrd_calc = XRDCalculator(debye_waller_factors={"Cu": 0.62, "Au": 0.58})
        xrd = xrd_calc.get_pattern(struct, scaled=False, two_theta_range=(0, 90))

        assert xrd.x == approx([42.58654814199049, 49.583339742779344, 72.7415377864004, 88.11241551060378])
        assert xrd.y == approx([2823918.0596514726, 1327811.2328217993, 770012.3450414365, 910580.4261507628])
        assert xrd.d_hkls == approx([2.1229169398102536, 1.8385, 1.3000158172114675, 1.1086572140124369])
        assert xrd.hkls == [
            [{"hkl": (1, 1, 1), "multiplicity": 8}],
            [{"hkl": (2, 0, 0), "multiplicity": 6}],
            [{"hkl": (2, 2, 0), "multiplicity": 12}],
            [{"hkl": (3, 1, 1), "multiplicity": 24}],
        ]

    def test_get_pattern_non_centrosymmetric(self):
        """Wurtzite ZnO (P6_3mc, non-centrosymmetric): F(hkl) is genuinely
        complex, so Friedel's law I(g) = I(-g) rests on F(-g) = F*(g) rather
        than the centrosymmetric special case of real F. Reference values were
        generated with the full-sphere (pre Friedel-halving) implementation
        and must not change. Multiplicities count both +l and -l reflections.
        """
        struct = Structure(
            Lattice.hexagonal(3.25, 5.207),
            ["Zn", "Zn", "O", "O"],
            [
                [1 / 3, 2 / 3, 0.0],
                [2 / 3, 1 / 3, 0.5],
                [1 / 3, 2 / 3, 0.3817],
                [2 / 3, 1 / 3, 0.8817],
            ],
        )
        xrd = XRDCalculator().get_pattern(struct, scaled=False, two_theta_range=(0, 90))

        assert len(xrd.x) == 13
        assert xrd.x[:4] == approx([31.793191209385288, 34.4481096267274, 36.28192404106642, 47.577360078761224])
        assert xrd.y[:4] == approx([138729.92514225902, 106194.79943567653, 270038.0917113368, 61581.215456583785])
        assert xrd.d_hkls[:4] == approx([2.814582562299426, 2.6035, 2.4760090064581086, 1.9112241235795175])
        assert xrd.hkls[:4] == [
            [{"hkl": (1, 0, -1, 0), "multiplicity": 6}],
            [{"hkl": (0, 0, 0, 2), "multiplicity": 2}],
            [{"hkl": (1, 0, -1, 1), "multiplicity": 12}],
            [{"hkl": (1, 0, -1, 2), "multiplicity": 12}],
        ]

    def test_get_pattern_chunked_matches_unchunked(self, monkeypatch):
        """Structure factors are accumulated in memory-bounded chunks of hkl
        rows; the result must not depend on the chunk size. Force one hkl row
        per chunk and compare against the effectively unchunked default.
        """
        struct = self.get_structure("LiFePO4")
        xrd_calc = XRDCalculator()
        ref = xrd_calc.get_pattern(struct, scaled=False, two_theta_range=(0, 90))

        monkeypatch.setattr(XRDCalculator, "PHASE_CHUNK_ENTRIES", 1)
        chunked = xrd_calc.get_pattern(struct, scaled=False, two_theta_range=(0, 90))

        assert chunked.x == approx(ref.x)
        assert chunked.y == approx(ref.y, rel=1e-12)
        assert chunked.d_hkls == approx(ref.d_hkls)
        assert chunked.hkls == ref.hkls

    def test_get_pattern_merge_regression(self):
        """
        Regression test for peak-merging behaviour. Assertion values were generated
        from the known-good reference implementation and must not change.
        """
        xrd_calc = XRDCalculator()

        # --- LiFePO4 ---
        struct = self.get_structure("LiFePO4")
        xrd = xrd_calc.get_pattern(struct, two_theta_range=(0, 90))

        assert len(xrd.hkls) == 434  # should have 434 peaks if no regression

        assert xrd.hkls[0:10] == [
            [{"hkl": (0, 1, 0), "multiplicity": 2}],
            [{"hkl": (0, 0, 2), "multiplicity": 2}],
            [{"hkl": (1, 0, -1), "multiplicity": 2}],
            [{"hkl": (1, 0, 1), "multiplicity": 2}],
            [{"hkl": (0, 1, -2), "multiplicity": 2}],
            [{"hkl": (0, 1, 2), "multiplicity": 2}],
            [{"hkl": (1, -1, 0), "multiplicity": 2}],
            [{"hkl": (1, 1, 0), "multiplicity": 2}],
            [{"hkl": (1, -1, 1), "multiplicity": 2}],
            [{"hkl": (1, 1, -1), "multiplicity": 2}],
        ]  # first 10 peaks are in this multiplicity.

        # --- Li10GeP2S12 ---
        struct = self.get_structure("Li10GeP2S12")
        xrd = xrd_calc.get_pattern(struct, two_theta_range=(0, 90))

        assert len(xrd.hkls) == 213
        assert xrd.hkls[0:10] == [
            [{"hkl": (1, 0, 1), "multiplicity": 8}],
            [{"hkl": (0, 0, 2), "multiplicity": 2}],
            [{"hkl": (1, 1, 0), "multiplicity": 4}],
            [{"hkl": (1, 0, 2), "multiplicity": 8}],
            [{"hkl": (1, 1, 2), "multiplicity": 8}],
            [{"hkl": (2, 0, 0), "multiplicity": 4}],
            [{"hkl": (2, 0, 1), "multiplicity": 8}],
            [{"hkl": (1, 0, 3), "multiplicity": 8}],
            [{"hkl": (2, 1, 1), "multiplicity": 16}],
            [{"hkl": (2, 0, 2), "multiplicity": 8}],
        ]
