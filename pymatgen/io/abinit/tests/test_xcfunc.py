# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.
from __future__ import unicode_literals, division, print_function

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinit.xcfunc import XcFunc

class LibxcFuncTest(PymatgenTest):

    def test_xcfunc_api(self):
        """Testing XcFunc API."""
        # GGA-PBE from ixc == 11
        ixc11 = XcFunc.from_abinit_ixc("11")
        assert ixc11.name == "PBE" and ixc11.type == "GGA"
        assert ixc11.name in XcFunc.known_names()

        d = {ixc11: ixc11.name}
        print(d)
        assert "PBE" in d
        assert ixc11 in d

        # Test if object can be serialized with Pickle.
        self.serialize_with_pickle(ixc11, test_eq=True)

        # Test if object supports MSONable
        # TODO
        #print("in test", type(ixc11.x), type(ixc11.c), type(ixc11.xc))
        #ixc11.x.as_dict()
        #self.assertMSONable(ixc11)

        # GGA-PBE from ixc given in abinit-libxc mode
        ixc_101130 = XcFunc.from_abinit_ixc("-101130")
        assert ixc_101130.name == "PBE" and ixc_101130.type == "GGA"
        assert ixc_101130 == ixc11

        # GGA-PBE built from name 
        gga_pbe = XcFunc.from_name("PBE")
        assert gga_pbe.name == "PBE" and gga_pbe.type == "GGA"
        assert ixc11 == gga_pbe

        #gga_pbe_libxc = XcFunc.from_type_name("GGA", "GGA_X_PBE+GGA_C_PBE")
        #assert ixc11 == gga_pbe_libxc
