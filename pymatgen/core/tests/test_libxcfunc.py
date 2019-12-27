# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.libxcfunc import LibxcFunc


class LibxcFuncTest(PymatgenTest):
    def test_libxcfunc_api(self):
        """Testing libxcfunc_api."""

        # LDA correlation: Hedin & Lundqvist
        xc = LibxcFunc.LDA_C_HL
        print(xc)
        assert not xc.is_x_kind and xc.is_c_kind and not xc.is_xc_kind
        assert xc.is_lda_family and not xc.is_gga_family
        print(xc.info_dict)

        assert xc.family in LibxcFunc.all_families()
        assert xc.kind in LibxcFunc.all_kinds()

        # Test if object can be serialized with Pickle.
        self.serialize_with_pickle(xc, test_eq=True)

        # Test if object supports MSONable
        self.assertMSONable(xc, test_if_subclass=False)
