from __future__ import annotations

from decimal import Decimal

import numpy as np
import pytest
from pytest import approx

from pymatgen.analysis.diffraction.core import DiffractionPattern
from pymatgen.analysis.diffraction.neutron import NDCalculator
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.phonon.thermal_displacements import ThermalDisplacementMatrices
from pymatgen.util.testing import MatSciTest

"""
These calculated values were verified with VESTA and FullProf.
"""

__author__ = "Yuta Suzuki"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Yuta Suzuki"
__email__ = "resnant@outlook.jp"
__date__ = "4/19/18"


class TestNDCalculator(MatSciTest):
    @pytest.mark.parametrize("two_theta_range", [(0, 5), (1, 5)])
    @pytest.mark.parametrize("scaled", [False, True])
    def test_get_pattern_empty_range(self, two_theta_range, scaled):
        # The first reflection is (100), at 2 theta = 25.45 degrees.
        structure = Structure(Lattice.cubic(3.5), ["Si"], [[0, 0, 0]])
        pattern = NDCalculator().get_pattern(structure, two_theta_range=two_theta_range, scaled=scaled)
        assert isinstance(pattern, DiffractionPattern)
        assert pattern.x.size == pattern.y.size == 0
        assert pattern.hkls == []
        assert pattern.d_hkls == []

    @pytest.mark.parametrize("two_theta_range", [(0, 5), (1, 5)])
    def test_empty_range_missing_scattering_coefficients(self, two_theta_range):
        structure = Structure(Lattice.cubic(3.5), ["Cm"], [[0, 0, 0]])
        with pytest.raises(ValueError, match="no scattering coefficients for Cm"):
            NDCalculator().get_pattern(structure, two_theta_range=two_theta_range)

    @pytest.mark.parametrize(
        "two_theta_range",
        [
            (5, 1),
            (-1, 5),
            (0, 181),
            (0, float("nan")),
            (0, float("inf")),
            (),
            (0,),
            (0, 5, 10),
            5,
            "05",
            ("0", 5),
            (None, 5),
            (0, 5j),
        ],
    )
    def test_get_pattern_invalid_range(self, two_theta_range):
        structure = Structure(Lattice.cubic(3.5), ["Si"], [[0, 0, 0]])
        with pytest.raises(
            ValueError, match="^two_theta_range must contain two finite real numbers satisfying 0 <= min <= max <= 180$"
        ):
            NDCalculator().get_pattern(structure, two_theta_range=two_theta_range)

    @pytest.mark.parametrize("two_theta_range", [None, (0, 90), [0, 90], np.array([0.0, 90.0])])
    def test_get_pattern_valid_range(self, two_theta_range):
        structure = Structure(Lattice.cubic(3.5), ["Si"], [[0, 0, 0]])
        pattern = NDCalculator().get_pattern(structure, two_theta_range=two_theta_range)
        assert len(pattern.x) > 0
        if two_theta_range is not None:
            assert len(pattern.x) == 9

    @pytest.mark.parametrize("range_factory", [lambda: (Decimal(0), Decimal(90)), lambda: iter([0, 90])])
    @pytest.mark.parametrize("scaled", [False, True])
    def test_get_pattern_compatible_range_types(self, range_factory, scaled):
        structure = Structure(Lattice.cubic(3.5), ["Si"], [[0, 0, 0]])
        calculator = NDCalculator()
        expected = calculator.get_pattern(structure, two_theta_range=(0, 90), scaled=scaled)
        pattern = calculator.get_pattern(structure, two_theta_range=range_factory(), scaled=scaled)
        np.testing.assert_array_equal(pattern.x, expected.x)
        np.testing.assert_array_equal(pattern.y, expected.y)
        assert pattern.hkls == expected.hkls
        assert pattern.d_hkls == expected.d_hkls

    @pytest.mark.parametrize("wavelength", [0, -1, float("inf"), float("nan")])
    def test_empty_range_invalid_wavelength(self, wavelength):
        structure = Structure(Lattice.cubic(3.5), ["Si"], [[0, 0, 0]])
        with pytest.raises(ValueError, match="wavelength must be finite and positive"):
            NDCalculator(wavelength=wavelength).get_pattern(structure, two_theta_range=(0, 5))

    def test_empty_structure(self):
        with pytest.raises(ValueError, match="structure must contain at least one site"):
            NDCalculator().get_pattern(Structure(Lattice.cubic(3.5), [], []), two_theta_range=(0, 5))

    @pytest.mark.parametrize("two_theta_range", [(0, 5), (1, 5)])
    def test_get_plot_empty_range(self, two_theta_range):
        import matplotlib.pyplot as plt

        structure = Structure(Lattice.cubic(3.5), ["Si"], [[0, 0, 0]])
        fig, ax = plt.subplots()
        try:
            assert NDCalculator().get_plot(structure, two_theta_range=two_theta_range, ax=ax) is ax
            assert len(ax.lines) == 0
            assert ax.get_xlabel()
            assert ax.get_ylabel()
        finally:
            plt.close(fig)

    @pytest.mark.parametrize("two_theta_range", [(179, 180), (0, 180), None])
    @pytest.mark.parametrize("scaled", [False, True])
    def test_get_pattern_backscattering_boundary(self, two_theta_range, scaled):
        # With a = 1 and wavelength = 2, the first shell is exactly at 180 degrees.
        structure = Structure(Lattice.cubic(1), ["Si"], [[0, 0, 0]])
        pattern = NDCalculator(wavelength=2).get_pattern(structure, two_theta_range=two_theta_range, scaled=scaled)
        assert pattern.x.size == pattern.y.size == 0
        assert pattern.hkls == []
        assert pattern.d_hkls == []

    @pytest.mark.parametrize("two_theta_range", [(0, 180), None])
    @pytest.mark.parametrize("scaled", [False, True])
    def test_backscattering_boundary_preserves_other_peaks(self, two_theta_range, scaled):
        structure = Structure(Lattice.tetragonal(1, 2), ["Si"], [[0, 0, 0]])
        calculator = NDCalculator(wavelength=2)
        expected = calculator.get_pattern(structure, two_theta_range=(0, 179), scaled=scaled)
        pattern = calculator.get_pattern(structure, two_theta_range=two_theta_range, scaled=scaled)
        assert len(expected.x) == 1
        np.testing.assert_array_equal(pattern.x, expected.x)
        np.testing.assert_array_equal(pattern.y, expected.y)
        assert pattern.hkls == expected.hkls
        assert pattern.d_hkls == expected.d_hkls

    @pytest.mark.parametrize("scaled", [False, True])
    def test_near_backscattering_peak_retained(self, scaled):
        structure = Structure(Lattice.cubic(1), ["Si"], [[0, 0, 0]])
        pattern = NDCalculator(wavelength=1.99999).get_pattern(structure, two_theta_range=(179, 180), scaled=scaled)
        assert len(pattern.x) == 1
        assert 179 < pattern.x[0] < 180
        assert np.isfinite(pattern.y[0])
        assert pattern.y[0] > 0

    def test_get_pattern(self):
        struct = self.get_structure("CsCl")
        c = NDCalculator(wavelength=1.54184)  # CuKa radiation
        pattern = c.get_pattern(struct, two_theta_range=(0, 90))
        # Check the first two peaks
        assert pattern.x[0] == approx(21.107738329639844)
        assert pattern.hkls[0] == [{"hkl": (1, 0, 0), "multiplicity": 6}]
        assert pattern.d_hkls[0] == approx(4.2089999999999996)
        assert pattern.x[1] == approx(30.024695921112777)
        assert pattern.hkls[1] == [{"hkl": (1, 1, 0), "multiplicity": 12}]
        assert pattern.d_hkls[1] == approx(2.976212442014178)

        struct = self.get_structure("LiFePO4")
        pattern = c.get_pattern(struct, two_theta_range=(0, 90))
        assert pattern.x[1] == approx(17.03504233621785)
        assert pattern.y[1] == approx(46.2985965)

        struct = self.get_structure("Li10GeP2S12")
        pattern = c.get_pattern(struct, two_theta_range=(0, 90))
        assert pattern.x[1] == approx(14.058274883353876)
        assert pattern.y[1] == approx(3.60588013)

        # Test a hexagonal structure.
        struct = self.get_structure("Graphite")
        pattern = c.get_pattern(struct, two_theta_range=(0, 90))

        ## REGRESSION TEST: The following values were verified
        #  with previous behaviour and should not change.
        assert pattern.x == approx(
            [
                26.21057350859598,
                42.249215447463456,
                44.39599754062177,
                50.40147067300912,
                53.934657179755774,
                59.41105634757505,
                70.7958623079758,
                77.25090657312492,
                83.23145147001156,
                84.40000457381223,
                85.72222578052532,
            ]
        )

        assert pattern.y == approx(
            [
                100.0,
                7.749071690702165,
                42.62382267,
                11.446571148133131,
                27.319449618794437,
                26.410629403328624,
                6.865301041009988,
                49.345768779989406,
                91.10964083311468,
                16.85276599985274,
                14.763509937143507,
            ]
        )
        assert len(pattern.hkls[0][0]) == approx(2)

        # Test an exception in case of the input element is
        # not in scattering length table.
        # This curium structure is just for test, not the actual structure.
        something = Structure(Lattice.cubic(a=1), ["Cm"], [[0, 0, 0]])
        with pytest.raises(
            ValueError,
            match=(
                "Unable to calculate ND pattern as "  # codespell:ignore ND
                "there is no scattering coefficients for Cm."
            ),
        ):
            pattern = c.get_pattern(something, two_theta_range=(0, 90))

        # Test with Debye-Waller factor
        struct = self.get_structure("Graphite")
        c = NDCalculator(wavelength=1.54184, debye_waller_factors={"C": 1})
        pattern = c.get_pattern(struct, two_theta_range=(0, 90))
        assert pattern.x[0] == approx(26.21057350859598)
        assert pattern.y[0] == approx(100)
        assert pattern.x[2] == approx(44.39599754)
        assert pattern.y[2] == approx(39.471514740)

    def test_get_pattern_non_centrosymmetric(self):
        """Tetragonal BaTiO3 (P4mm, non-centrosymmetric) with Ti, whose
        neutron scattering length is negative (-3.438 fm): F(hkl) is genuinely
        complex, so Friedel's law I(g) = I(-g) rests on F(-g) = F*(g) rather
        than the centrosymmetric special case of real F. Reference values were
        generated with the full-sphere (pre Friedel-halving) implementation
        and must not change. Multiplicities count both +l and -l reflections,
        e.g. (001) has multiplicity 2 although P4mm does not relate +l to -l.
        """
        struct = Structure(
            Lattice.tetragonal(3.9945, 4.0334),
            ["Ba", "Ti", "O", "O", "O"],
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.4820],
                [0.5, 0.5, 0.0160],
                [0.5, 0.0, 0.5150],
                [0.0, 0.5, 0.5150],
            ],
        )
        calc = NDCalculator(wavelength=1.54184)
        pattern = calc.get_pattern(struct, scaled=False, two_theta_range=(0, 90))

        assert len(pattern.x) == 27
        assert pattern.x[:4] == approx([22.037944670652465, 22.255284684116262, 31.52181770727521, 31.678181395101802])
        assert pattern.y[:4] == approx([453.13576368164007, 800.8355861573295, 1916.4711408968653, 971.0175541753043])
        assert pattern.d_hkls[:4] == approx([4.0334, 3.9945, 2.838191301738583, 2.824538037449664])
        assert pattern.hkls[:4] == [
            [{"hkl": (0, 0, 1), "multiplicity": 2}],
            [{"hkl": (1, 0, 0), "multiplicity": 4}],
            [{"hkl": (1, 0, 1), "multiplicity": 8}],
            [{"hkl": (1, 1, 0), "multiplicity": 4}],
        ]

    def test_get_pattern_chunked_matches_unchunked(self, monkeypatch):
        """Structure factors are accumulated in memory-bounded chunks of hkl
        rows; the result must not depend on the chunk size. Force one hkl row
        per chunk and compare against the effectively unchunked default.
        """
        struct = self.get_structure("LiFePO4")
        calc = NDCalculator(wavelength=1.54184, debye_waller_factors={"Fe": 0.4})
        ref = calc.get_pattern(struct, scaled=False, two_theta_range=(0, 90))

        monkeypatch.setattr(NDCalculator, "PHASE_CHUNK_ENTRIES", 1)
        chunked = calc.get_pattern(struct, scaled=False, two_theta_range=(0, 90))

        assert chunked.x == approx(ref.x)
        assert chunked.y == approx(ref.y, rel=1e-12)
        assert chunked.d_hkls == approx(ref.d_hkls)
        assert chunked.hkls == ref.hkls

    def test_get_plot(self):
        struct = self.get_structure("Graphite")
        c = NDCalculator(wavelength=1.54184, debye_waller_factors={"C": 1})
        c.get_plot(struct, two_theta_range=(0, 90))

    def test_get_pattern_isotropic_thermal_displacements(self, monkeypatch):
        """An isotropic U gives the same pattern as debye_waller_factors with B = 8 pi^2 U."""
        struct = self.get_structure("LiFePO4")
        struct.replace_species({"Li": {"Li": 0.5, "Na": 0.5}})
        u_iso = {"Li": 0.02, "Na": 0.02, "Fe": 0.005, "P": 0.004, "O": 0.008}
        b_factors = {el: 8 * np.pi**2 * u for el, u in u_iso.items()}
        ref = NDCalculator(debye_waller_factors=b_factors).get_pattern(struct, scaled=False)
        u_red = ThermalDisplacementMatrices.get_reduced_matrix(
            [u_iso[site.species.elements[0].symbol] * np.eye(3) for site in struct]
        )
        struct = ThermalDisplacementMatrices(u_red, struct, temperature=None).to_structure_with_site_properties_Ucif()
        assert NDCalculator().get_pattern(struct, scaled=False).y == approx(ref.y)
        # one hkl row per chunk, so the factors are also checked across chunks
        monkeypatch.setattr(NDCalculator, "PHASE_CHUNK_ENTRIES", 1)
        pattern = NDCalculator().get_pattern(struct, scaled=False)
        assert pattern.x == approx(ref.x)
        assert pattern.y == approx(ref.y)

    @pytest.mark.parametrize("kwargs", [{"symprec": 0.1}, {"debye_waller_factors": {"Cs": 0.1}}])
    def test_get_pattern_thermal_displacements_conflicts(self, kwargs):
        struct = self.get_structure("CsCl")
        struct.add_site_property("magmom", [0, 0])
        NDCalculator(**kwargs).get_pattern(struct)
        u_red = ThermalDisplacementMatrices.get_reduced_matrix([0.01 * np.eye(3)] * len(struct))
        struct = ThermalDisplacementMatrices(u_red, struct, temperature=None).to_structure_with_site_properties_Ucif()
        with pytest.raises(ValueError, match="Neither symprec nor debye_waller_factors can be used with the U11_cif"):
            NDCalculator(**kwargs).get_pattern(struct)
