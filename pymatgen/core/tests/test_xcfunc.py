# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from __future__ import unicode_literals, division, print_function

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.xcfunc import XcFunc


class LibxcFuncTest(PymatgenTest):

    def test_xcfunc_api(self):
        """Testing XcFunc API."""
        # Aliases should be unique
        assert len(XcFunc.aliases()) == len(set(XcFunc.aliases()))

        # LDA-Teter
        ixc_1 = XcFunc.from_abinit_ixc(1)
        print(ixc_1)
        assert ixc_1.type == "LDA"
        assert ixc_1.name == "LDA_XC_TETER93"
        assert ixc_1 == ixc_1
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

        # GGA-PBE from ixc == 11 (in aliases)
        ixc_11 = XcFunc.from_abinit_ixc(11)
        assert ixc_11.type == "GGA" and ixc_11.name == "PBE"
        assert ixc_11.name in XcFunc.aliases()
        assert ixc_1 != ixc_11
        # Test asxc
        assert XcFunc.asxc(ixc_11) is ixc_11
        assert XcFunc.asxc("PBE") == ixc_11

        d = {ixc_11: ixc_11.name}
        print(d)
        assert "PBE" in d
        assert ixc_11 in d

        # Test if object can be serialized with Pickle.
        self.serialize_with_pickle(ixc_11, test_eq=True)

        # Test if object supports MSONable
        # TODO
        #print("in test", type(ixc_11.x), type(ixc_11.c), type(ixc_11.xc))
        #ixc_11.x.as_dict()
        #self.assertMSONable(ixc_11)

        # GGA-PBE from ixc given in abinit-libxc mode
        ixc_101130 = XcFunc.from_abinit_ixc(-101130)
        assert ixc_101130.type == "GGA" and ixc_101130.name == "PBE"
        assert ixc_101130 == ixc_11

        # GGA-PBE built from name
        gga_pbe = XcFunc.from_name("PBE")
        assert gga_pbe.type == "GGA" and gga_pbe.name == "PBE"
        assert ixc_11 == gga_pbe

        # Use X from GGA and C from LDA!
        unknown_xc = XcFunc.from_name("GGA_X_PBE+ LDA_C_PW")
        assert unknown_xc not in XcFunc.aliases()
        assert unknown_xc.type == "GGA+LDA"
        assert unknown_xc.name == "GGA_X_PBE+LDA_C_PW"

        gga_pbe = XcFunc.from_type_name("GGA", "GGA_X_PBE+GGA_C_PBE")
        assert gga_pbe.type == "GGA" and gga_pbe.name == "PBE"
        assert str(gga_pbe) == "PBE"


