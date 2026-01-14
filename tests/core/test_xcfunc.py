from __future__ import annotations

import pytest

from pymatgen.core.xcfunc import XcFunc
from pymatgen.util.testing import MatSciTest


class TestLibxcFunc(MatSciTest):
    def setup_method(self) -> None:
        self.ixc_11 = XcFunc.from_abinit_ixc(11)

    def test_aliases(self):
        # Aliases should be unique
        assert len(XcFunc.aliases()) == len(set(XcFunc.aliases()))

    def test_lda(self):
        # LDA-Teter
        ixc_1 = XcFunc.from_abinit_ixc(1)
        assert ixc_1.type == "LDA"
        assert ixc_1.name == "LDA_XC_TETER93"
        assert ixc_1 == ixc_1  # test __eq__  # noqa: PLR0124
        assert ixc_1 == "LDA_XC_TETER93"
        assert ixc_1 != "PBE"
        assert ixc_1.name not in XcFunc.aliases()
        assert ixc_1 == XcFunc.from_name(ixc_1.name)

        # LDA-PW (in aliases)
        ixc_7 = XcFunc.from_abinit_ixc(7)
        assert ixc_7.type == "LDA"
        assert ixc_7.name == "PW"
        assert ixc_7.name in XcFunc.aliases()
        assert ixc_7.name == XcFunc.from_name(ixc_7.name)
        assert ixc_7 != ixc_1

    def test_gga_pbe(self):
        # GGA-PBE from ixc == 11 (in aliases)
        ixc_1 = XcFunc.from_abinit_ixc(1)
        ixc_11 = XcFunc.from_abinit_ixc(11)
        assert ixc_11.type == "GGA"
        assert ixc_11.name == "PBE"
        assert ixc_11.name in XcFunc.aliases()
        assert ixc_1 != ixc_11
        # Test asxc
        assert XcFunc.asxc(ixc_11) is ixc_11
        assert XcFunc.asxc("PBE") == ixc_11

        dct = {ixc_11: ixc_11.name}
        assert "PBE" in dct
        assert ixc_11 in dct

    def test_pickle_serialize(self):
        # Test if object can be serialized with Pickle
        self.serialize_with_pickle(self.ixc_11)

    @pytest.mark.xfail(reason="TODO:")
    def test_msonable(self):
        # Test if object supports MSONable
        self.ixc_11.x.as_dict()
        self.assert_msonable(self.ixc_11)

    def test_from(self):
        # GGA-PBE from ixc given in abinit-libxc mode
        ixc_101130 = XcFunc.from_abinit_ixc(-101130)
        assert ixc_101130.type == "GGA"
        assert ixc_101130.name == "PBE"
        assert ixc_101130 == self.ixc_11

        # GGA-PBE built from name
        gga_pbe = XcFunc.from_name("PBE")
        assert gga_pbe.type == "GGA"
        assert gga_pbe.name == "PBE"
        assert self.ixc_11 == gga_pbe

        # Use X from GGA and C from LDA!
        unknown_xc = XcFunc.from_name("GGA_X_PBE+ LDA_C_PW")
        assert unknown_xc not in XcFunc.aliases()
        assert unknown_xc.type == "GGA+LDA"
        assert unknown_xc.name == "GGA_X_PBE+LDA_C_PW"

        gga_pbe = XcFunc.from_type_name("GGA", "GGA_X_PBE+GGA_C_PBE")
        assert gga_pbe.type == "GGA"
        assert gga_pbe.name == "PBE"
        assert str(gga_pbe) == "PBE"
