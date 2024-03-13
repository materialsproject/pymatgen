from __future__ import annotations

from pymatgen.core.libxcfunc import LibxcFunc
from pymatgen.util.testing import PymatgenTest


class TestLibxcFunc(PymatgenTest):
    def test_libxcfunc_api(self):
        """Testing libxcfunc_api."""
        # LDA correlation: Hedin & Lundqvist
        xc = LibxcFunc.LDA_C_HL
        assert not xc.is_x_kind
        assert xc.is_c_kind
        assert not xc.is_xc_kind
        assert xc.is_lda_family
        assert not xc.is_gga_family

        assert xc.family in LibxcFunc.all_families()
        assert xc.kind in LibxcFunc.all_kinds()

        # Test if object can be serialized with Pickle.
        self.serialize_with_pickle(xc)

        # Test if object supports MSONable
        self.assert_msonable(xc, test_is_subclass=False)

    def test_repr(self):
        """Test LibxcFunc.__repr__"""
        xc = LibxcFunc.LDA_C_HL
        assert repr(xc) == "LibxcFunc(name='LDA_C_HL', kind='CORRELATION', family='LDA')"
