from __future__ import annotations

import unittest

import numpy as np
import pytest
from pytest import approx

from pymatgen.core.lattice import Lattice
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.groups import SYMM_DATA, PointGroup, SpaceGroup

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Virtual Lab"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "ongsp@ucsd.edu"
__date__ = "4/10/14"

ORDERED_SYMBOLS = (
    "P1 P-1 P121 P12_11 C121 P1m1 P1c1 C1m1 C1c1 P12/m1 P12_1/m1 C12/m1 P12/c1 P12_1/c1 C12/c1 P222 P222_1"
    " P2_12_12 P2_12_121 C222_1 C222 F222 I222 I2_12_121 Pmm2 Pmc2_1 Pcc2 Pma2 Pca2_1 Pnc2 Pmn2_1 Pba2 Pna2_1 Pnn2 "
    "Cmm2 Cmc2_1 Ccc2 Amm2 Aem2 Ama2 Aea2 Fmm2 Fdd2 Imm2 Iba2 Ima2 Pmmm Pnnn1 Pccm Pban1 Pmma Pnna Pmna Pcca Pbam "
    "Pccn Pbcm Pnnm Pmmn1 Pbcn Pbca Pnma Cmcm Cmce Cmmm Cccm Cmme Ccce1 Fmmm Fddd1 Immm Ibam Ibca Imma P4 P4_1 P4_2 "
    "P4_3 I4 I4_1 P-4 I-4 P4/m P4_2/m P4/n1 P4_2/n I4/m I4_1/a P422 P42_12 P4_122 P4_12_12 P4_222 P4_22_12 P4_322 "
    "P4_32_12 I422 I4_122 P4mm P4bm P4_2cm P4_2nm P4cc P4nc P4_2mc P4_2bc I4mm I4cm I4_1md I4_1cd P-42m P-42c P-42_1m "
    "P-42_1c P-4m2 P-4c2 P-4b2 P-4n2 I-4m2 I-4c2 I-42m I-42d P4/mmm P4/mcc P4/nbm1 P4/nnc1 P4/mbm P4/mnc P4/nmm1 "
    "P4/ncc1 P4_2/mmc P4_2/mcm P4_2/nbc P4_2/nnm P4_2/mbc P4_2/mnm P4_2/nmc P4_2/ncm I4/mmm I4/mcm I4_1/amd I4_1/acd "
    "P3 P3_1 P3_2 R3H P-3 R-3H P312 P321 P3_112 P3_121 P3_212 P3_221 R32H P3m1 P31m P3c1 P31c R3mH R3cH P-31m P-31c "
    "P-3m1 P-3c1 R-3mH R-3cH P6 P6_1 P6_5 P6_2 P6_4 P6_3 P-6 P6/m P6_3/m P622 P6_122 P6_522 P6_222 P6_422 P6_322 "
    "P6mm P6cc P6_3cm P6_3mc P-6m2 P-6c2 P-62m P-62c P6/mmm P6/mcc P6_3/mcm P6_3/mmc P23 F23 I23 P2_13 I2_13 Pm-3 "
    "Pn-31 Fm-3 Fd-31 Im-3 Pa-3 Ia-3 P432 P4_232 F432 F4_132 I432 P4_332 P4_132 I4_132 P-43m F-43m I-43m P-43n F-43c "
    "I-43d Pm-3m Pn-3n1 Pm-3n Pn-3m1 Fm-3m Fm-3c Fd-3m1 Fd-3c1 Im-3m Ia-3d"
).split()


class TestPointGroup(unittest.TestCase):
    def test_order(self):
        orders = {"mmm": 8, "432": 24, "-6m2": 12}
        for key, val in orders.items():
            assert len(PointGroup(key).symmetry_ops) == val

    def test_get_orbit(self):
        pg_mmm = PointGroup("mmm")
        assert len(pg_mmm.get_orbit([0.1, 0.1, 0.1])) == 8
        assert len(pg_mmm.get_orbit([0, 0, 0.1])) == 2
        assert len(pg_mmm.get_orbit([1.2, 1.2, 1])) == 8

    def test_is_sub_super_group(self):
        pg_mmm = PointGroup("mmm")
        pg_mm2 = PointGroup("mm2")
        pg_222 = PointGroup("222")
        pg_4 = PointGroup("4")
        assert pg_mmm.is_supergroup(pg_mm2)
        assert pg_mm2.is_subgroup(pg_mmm)
        assert pg_mmm.is_supergroup(pg_222)
        assert not pg_mmm.is_supergroup(pg_4)
        pg_m3m = PointGroup("m-3m")
        pg_6mmm = PointGroup("6/mmm")
        pg_3m = PointGroup("-3m")
        # TODO: Fix the test below.
        # assert pg3m.is_subgroup(pgm3m)
        assert pg_3m.is_subgroup(pg_6mmm)
        assert not pg_m3m.is_supergroup(pg_6mmm)


class TestSpaceGroup(unittest.TestCase):
    def test_renamed_e_symbols(self):
        assert SpaceGroup.from_int_number(64).symbol == "Cmce"

        for sym, num in (("Aem2", 39), ("Aea2", 41), ("Cmce", 64), ("Cmme", 67), ("Ccce", 68)):
            assert SpaceGroup(sym).int_number == num

    def test_abbrev_symbols(self):
        sg = SpaceGroup("P2/c")
        assert sg.int_number == 13
        sg = SpaceGroup("R-3mH")
        assert sg.int_number == 166

    def test_attr(self):
        sg = SpaceGroup("Fm-3m")
        assert sg.full_symbol == "F4/m-32/m"
        assert sg.point_group == "m-3m"

    def test_point_group_is_set(self):
        for num in range(1, 231):
            sg = SpaceGroup.from_int_number(num)
            assert hasattr(sg, "point_group")

        for symbol in SYMM_DATA["space_group_encoding"]:
            sg = SpaceGroup(symbol)
            assert hasattr(sg, "point_group")

    def test_full_symbols(self):
        sg = SpaceGroup("P2/m2/m2/m")
        assert sg.symbol == "Pmmm"

    def test_order_symm_ops(self):
        for name in SpaceGroup.SG_SYMBOLS:
            sg = SpaceGroup(name)
            assert len(sg.symmetry_ops) == sg.order

    def test_get_settings(self):
        assert SpaceGroup.get_settings("Fm-3m") == {"Fm-3m(a-1/4,b-1/4,c-1/4)", "Fm-3m"}
        assert SpaceGroup.get_settings("Pmmn") == {
            "Pmmn",
            "Pmnm:1",
            "Pnmm:2",
            "Pmnm:2",
            "Pnmm",
            "Pnmm:1",
            "Pmmn:1",
            "Pmnm",
            "Pmmn:2",
        }
        assert SpaceGroup.get_settings("Pmna") == {"Pnmb", "Pman", "Pncm", "Pmna", "Pcnm", "Pbmn"}

    def test_crystal_system(self):
        sg = SpaceGroup("R-3c")
        assert sg.crystal_system == "trigonal"
        sg = SpaceGroup("R-3cH")
        assert sg.crystal_system == "trigonal"

    def test_get_orbit(self):
        sg = SpaceGroup("Fm-3m")
        p = np.random.randint(0, 100 + 1, size=(3,)) / 100
        assert len(sg.get_orbit(p)) <= sg.order

    def test_get_orbit_and_generators(self):
        sg = SpaceGroup("Fm-3m")
        p = np.random.randint(0, 100 + 1, size=(3,)) / 100
        orbit, generators = sg.get_orbit_and_generators(p)
        assert len(orbit) <= sg.order
        pp = generators[0].operate(orbit[0])
        assert p[0] == approx(pp[0])
        assert p[1] == approx(pp[1])
        assert p[2] == approx(pp[2])

    def test_is_compatible(self):
        cubic = Lattice.cubic(1)
        hexagonal = Lattice.hexagonal(1, 2)
        rhom = Lattice.rhombohedral(3, 80)
        tet = Lattice.tetragonal(1, 2)
        ortho = Lattice.orthorhombic(1, 2, 3)
        sg = SpaceGroup("Fm-3m")
        assert sg.is_compatible(cubic)
        assert not sg.is_compatible(hexagonal)
        sg = SpaceGroup("R-3m:H")
        assert not sg.is_compatible(cubic)
        assert sg.is_compatible(hexagonal)
        sg = SpaceGroup("R-3m:R")
        assert sg.is_compatible(cubic)
        assert sg.is_compatible(rhom)
        assert not sg.is_compatible(hexagonal)
        sg = SpaceGroup("Pnma")
        assert sg.is_compatible(cubic)
        assert sg.is_compatible(tet)
        assert sg.is_compatible(ortho)
        assert not sg.is_compatible(rhom)
        assert not sg.is_compatible(hexagonal)
        sg = SpaceGroup("P12/c1")
        assert sg.is_compatible(cubic)
        assert sg.is_compatible(tet)
        assert sg.is_compatible(ortho)
        assert not sg.is_compatible(rhom)
        assert not sg.is_compatible(hexagonal)
        sg = SpaceGroup("P-1")
        assert sg.is_compatible(cubic)
        assert sg.is_compatible(tet)
        assert sg.is_compatible(ortho)
        assert sg.is_compatible(rhom)
        assert sg.is_compatible(hexagonal)
        sg = SpaceGroup("Pmmn:2")
        assert sg.is_compatible(cubic)
        assert sg.is_compatible(tet)
        assert sg.is_compatible(ortho)
        assert not sg.is_compatible(rhom)
        assert not sg.is_compatible(hexagonal)

        sg = SpaceGroup.from_int_number(165)
        assert not sg.is_compatible(cubic)
        assert not sg.is_compatible(tet)
        assert not sg.is_compatible(ortho)
        assert not sg.is_compatible(rhom)
        assert sg.is_compatible(hexagonal)

    def test_symmops(self):
        sg = SpaceGroup("Pnma")
        op = SymmOp.from_rotation_and_translation([[1, 0, 0], [0, -1, 0], [0, 0, -1]], [0.5, 0.5, 0.5])
        assert op in sg.symmetry_ops

    def test_other_settings(self):
        sg = SpaceGroup("Pbnm")
        assert sg.int_number == 62
        assert sg.order == 8
        with pytest.raises(ValueError, match="Bad international symbol 'hello'"):
            SpaceGroup("hello")

    def test_subgroup_supergroup(self):
        assert SpaceGroup("Pma2").is_subgroup(SpaceGroup("Pccm"))
        assert not SpaceGroup.from_int_number(229).is_subgroup(SpaceGroup.from_int_number(230))

    def test_hexagonal(self):
        for num in (146, 148, 155, 160, 161, 166, 167):
            sg = SpaceGroup.from_int_number(num, hexagonal=False)
            assert not sg.symbol.endswith("H")

    def test_string(self):
        sg = SpaceGroup("R-3c")
        assert sg.to_latex_string() == r"R$\overline{3}$cH"
        sg = SpaceGroup("P6/mmm")
        assert sg.to_latex_string() == "P6/mmm"
        sg = SpaceGroup("P4_1")
        assert sg.to_unicode_string() == "P4₁"

    def test_repr(self):
        for num in range(1, 231):
            sg = SpaceGroup.from_int_number(num)
            symbol = ORDERED_SYMBOLS[num - 1]
            assert repr(sg) in f"SpaceGroup({symbol=})"

    def test_raises_on_bad_int_number(self):
        for num in (-5, 0, 231, 1000):
            with pytest.raises(ValueError, match=f"International number must be between 1 and 230, got {num}"):
                SpaceGroup.from_int_number(num)
