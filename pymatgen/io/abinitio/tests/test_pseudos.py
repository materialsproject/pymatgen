# coding: utf-8

from __future__ import unicode_literals, division, print_function

"""
Created on Fri Mar  8 23:14:02 CET 2013
"""

import os.path
import collections

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinitio.pseudos import *

_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files', "abinitio")


def ref_file(filename):
    return os.path.join(_test_dir, filename)


def ref_files(*filenames):
    return list(map(ref_file, filenames))


class PseudoTestCase(PymatgenTest):

    def setUp(self):
        nc_pseudo_fnames = collections.defaultdict(list)
        nc_pseudo_fnames["Si"] = ref_files("14si.pspnc",  "14si.4.hgh", "14-Si.LDA.fhi")

        self.nc_pseudos = collections.defaultdict(list)

        for (symbol, fnames) in nc_pseudo_fnames.items():
            for fname in fnames:
                root, ext = os.path.splitext(fname)
                pseudo = Pseudo.from_file(ref_file(fname))
                self.nc_pseudos[symbol].append(pseudo)

                # Save the pseudo as instance attribute whose name 
                # is constructed with the rule: symbol_ppformat
                attr_name = symbol + "_" + ext[1:]
                if hasattr(self, attr_name):
                    raise RuntimError("self has already the attribute %s" % attr_name)

                setattr(self, attr_name, pseudo)

    def test_nc_pseudos(self):
        """Test norm-conserving pseudopotentials"""
        for (symbol, pseudos) in self.nc_pseudos.items():
            for pseudo in pseudos:
                print(repr(pseudo))
                print(pseudo)
                self.assertTrue(pseudo.isnc)
                self.assertFalse(pseudo.ispaw)
                self.assertEqual(pseudo.Z, 14)
                self.assertEqual(pseudo.symbol, symbol)
                self.assertEqual(pseudo.Z_val, 4)
                self.assertGreaterEqual(pseudo.nlcc_radius, 0.0)
                print(pseudo.as_dict())

                self.assertPMGSONable(pseudo)

                # Test pickle
                self.serialize_with_pickle(pseudo, test_eq=False)

        # HGH pseudos
        pseudo = self.Si_hgh
        self.assertFalse(pseudo.has_nlcc)
        self.assertEqual(pseudo.l_max, 1)
        self.assertEqual(pseudo.l_local, 0)

        assert self.Si_hgh == self.Si_hgh

        # TM pseudos
        pseudo = self.Si_pspnc
        self.assertTrue(pseudo.has_nlcc)
        self.assertEqual(pseudo.l_max, 2)
        self.assertEqual(pseudo.l_local, 2)

        assert self.Si_hgh != self.Si_pspnc

        # FHI pseudos
        pseudo = self.Si_fhi
        self.assertFalse(pseudo.has_nlcc)
        self.assertEqual(pseudo.l_max, 3)
        self.assertEqual(pseudo.l_local, 2)
        
        # Test PseudoTable.
        table = PseudoTable(self.nc_pseudos["Si"])
        print(repr(table))
        print(table)
        self.assertTrue(table.allnc)
        self.assertTrue(not table.allpaw)
        self.assertFalse(not table.is_complete)
        assert len(table) == 3
        assert len(table[14]) == 3
        assert len(table.select_symbols("Si")) == 3
        assert table.zlist == [14]

        # Test pickle
        self.serialize_with_pickle(table, test_eq=False)

    def test_pawxml_pseudos(self):
        """Test O.GGA_PBE-JTH-paw.xml."""
        oxygen = Pseudo.from_file(ref_file("O.GGA_PBE-JTH-paw.xml"))
        print(repr(oxygen))
        print(oxygen)
        print(oxygen.as_dict())

        self.assertTrue(oxygen.ispaw)
        self.assertTrue(oxygen.symbol == "O" and 
                       (oxygen.Z, oxygen.core, oxygen.valence) == (8, 2, 6),
                        oxygen.Z_val == 6,
                       )

        self.assert_almost_equal(oxygen.paw_radius, 1.4146523028)

        # Test pickle
        new_objs = self.serialize_with_pickle(oxygen, test_eq=False)

        for o in new_objs:
            print(repr(o))
            print(o)
                                                                                 
            self.assertTrue(o.ispaw)
            self.assertTrue(o.symbol == "O" and 
                           (o.Z, o.core, o.valence) == (8, 2, 6),
                            o.Z_val == 6,
                           )
                                                                                 
            self.assert_almost_equal(o.paw_radius, 1.4146523028)

    def test_oncvpsp_pseudo(self):
        """
        Test the ONCVPSP Ge pseudo
        """
        ger = Pseudo.from_file(ref_file("ge.oncvpsp"))
        print(repr(ger))
        print(ger)
        print(ger.as_dict())

        self.assertTrue(ger.symbol == "Ge")
        self.assert_equal(ger.Z, 32.0)
        self.assert_equal(ger.Z_val, 4.0)
        self.assertTrue(ger.isnc)
        self.assertFalse(ger.ispaw)
        self.assert_equal(ger.l_max, 2)
        self.assert_equal(ger.l_local, 4)
        self.assert_equal(ger.rcore, None)
        self.assertFalse(ger.has_dojo_report)

    def test_oncvpsp_dojo_report(self):
        """
        Test the dojo report
        """
        plot = True

        try:
            from matplotlib.figure import Figure as Fig
        except ImportError:
            Fig = None
            plot = False

        h_wdr = Pseudo.from_file(ref_file("H-wdr.oncvpsp"))

        print(repr(h_wdr))
        print(h_wdr.as_dict())

        self.assertTrue(h_wdr.symbol == "H")
        self.assertTrue(h_wdr.has_dojo_report)

        report = h_wdr.read_dojo_report()
        self.assert_equal(report.check(), {})

        if plot:
            self.assertIsInstance(report.plot_deltafactor_convergence(show=False), Fig)
            self.assertIsInstance(report.plot_deltafactor_eos(show=False), Fig)
            self.assertIsInstance(report.plot_etotal_vs_ecut(show=False), Fig)
            self.assertIsInstance(report.plot_gbrv_convergence(show=False), Fig)
            self.assertIsInstance(report.plot_gbrv_eos('bcc', show=False), Fig)
            self.assertIsInstance(report.plot_gbrv_eos('fcc', show=False), Fig)
            self.assertIsInstance(report.plot_phonon_convergence(show=False), Fig)

       # self.assertFalse(report.has_exceptions())


class PseudoTableTest(PymatgenTest):

    def test_methods(self):
        """Test PseudoTable methods"""
        table = PseudoTable(ref_files("14si.pspnc",  "14si.4.hgh", "14-Si.LDA.fhi"))
        print(table)
        assert len(table) == 3
        for pseudo in table:
            assert pseudo.isnc
        assert table.allnc and not table.allpaw
        assert table.zlist == [14]

        # Data persistence
        self.serialize_with_pickle(table, test_eq=False)

        #d = table.as_dict()
        #PseudoTable.from_dict(d)
        #self.assertPMGSONable(table)

        selected = table.select_symbols("Si")
        assert len(selected) == len(table) and selected.__class__ is table.__class__

        with self.assertRaises(ValueError):
            table.pseudos_with_symbols("Si")


if __name__ == "__main__":
    import unittest
    unittest.main()
