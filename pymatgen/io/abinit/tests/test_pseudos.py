# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import unicode_literals, division, print_function

import os.path
import collections
import numpy as np

from pymatgen.util.testing import PymatgenTest
from pymatgen.io.abinit.pseudos import *


_test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files', "abinit")

def ref_file(filename):
    return os.path.join(_test_dir, filename)


def ref_files(*filenames):
    return list(map(ref_file, filenames))


class PseudoTestCase(PymatgenTest):

    def setUp(self):
        nc_pseudo_fnames = collections.defaultdict(list)
        nc_pseudo_fnames["Si"] = ref_files("14si.pspnc",  "14si.4.hgh", "14-Si.LDA.fhi")

        self.nc_pseudos = collections.defaultdict(list)

        for symbol, fnames in nc_pseudo_fnames.items():
            for fname in fnames:
                root, ext = os.path.splitext(fname)
                pseudo = Pseudo.from_file(fname)
                self.nc_pseudos[symbol].append(pseudo)

                # Save the pseudo as instance attribute whose name
                # is constructed with the rule: symbol_ppformat
                attr_name = symbol + "_" + ext[1:]
                if hasattr(self, attr_name):
                    raise RuntimeError("self has already the attribute %s" % attr_name)

                setattr(self, attr_name, pseudo)

    def test_nc_pseudos(self):
        """Test norm-conserving pseudopotentials"""
        for symbol, pseudos in self.nc_pseudos.items():
            for pseudo in pseudos:
                print(repr(pseudo))
                print(pseudo)
                self.assertTrue(pseudo.isnc)
                self.assertFalse(pseudo.ispaw)
                self.assertEqual(pseudo.Z, 14)
                self.assertEqual(pseudo.symbol, symbol)
                self.assertEqual(pseudo.Z_val, 4)
                self.assertGreaterEqual(pseudo.nlcc_radius, 0.0)

                # Test pickle
                self.serialize_with_pickle(pseudo, test_eq=False)

                # Test MSONable
                #print(pseudo.as_dict())
                self.assertMSONable(pseudo)

        # HGH pseudos
        pseudo = self.Si_hgh
        self.assertFalse(pseudo.has_nlcc)
        self.assertEqual(pseudo.l_max, 1)
        self.assertEqual(pseudo.l_local, 0)
        assert not pseudo.supports_soc
        assert self.Si_hgh.md5 is not None
        assert self.Si_hgh == self.Si_hgh

        # TM pseudos
        pseudo = self.Si_pspnc
        self.assertTrue(pseudo.has_nlcc)
        self.assertEqual(pseudo.l_max, 2)
        self.assertEqual(pseudo.l_local, 2)
        assert not pseudo.supports_soc
        assert self.Si_hgh != self.Si_pspnc

        # FHI pseudos
        pseudo = self.Si_fhi
        self.assertFalse(pseudo.has_nlcc)
        self.assertEqual(pseudo.l_max, 3)
        self.assertEqual(pseudo.l_local, 2)
        assert not pseudo.supports_soc

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

        assert oxygen.xc.type == "GGA" and oxygen.xc.name == "PBE"
        assert oxygen.supports_soc
        assert oxygen.md5 is not None
        self.assert_almost_equal(oxygen.paw_radius, 1.4146523028)

        # Test pickle
        new_objs = self.serialize_with_pickle(oxygen, test_eq=False)
        # Test MSONable
        self.assertMSONable(oxygen)

        for o in new_objs:
            print(repr(o))
            print(o)

            self.assertTrue(o.ispaw)
            self.assertTrue(o.symbol == "O" and
                           (o.Z, o.core, o.valence) == (8, 2, 6),
                            o.Z_val == 6,
                           )

            self.assert_almost_equal(o.paw_radius, 1.4146523028)

    def test_oncvpsp_pseudo_sr(self):
        """
        Test the ONCVPSP Ge pseudo (scalar relativistic version).
        """
        ger = Pseudo.from_file(ref_file("ge.oncvpsp"))
        print(repr(ger))
        print(ger)
        print(ger.as_dict())
        ger.as_tmpfile()

        self.assertTrue(ger.symbol == "Ge")
        self.assert_equal(ger.Z, 32.0)
        self.assert_equal(ger.Z_val, 4.0)
        self.assertTrue(ger.isnc)
        self.assertFalse(ger.ispaw)
        self.assert_equal(ger.l_max, 2)
        self.assert_equal(ger.l_local, 4)
        self.assert_equal(ger.rcore, None)
        assert not ger.supports_soc

        # Data persistence
        self.serialize_with_pickle(ger, test_eq=False)
        self.assertMSONable(ger)

    def test_oncvpsp_pseudo_fr(self):
        """
        Test the ONCVPSP Pb pseudo (relativistic version with SO).
        """
        pb = Pseudo.from_file(ref_file("Pb-d-3_r.psp8"))
        print(repr(pb))
        print(pb)
        #print(pb.as_dict())
        #pb.as_tmpfile()

        # Data persistence
        self.serialize_with_pickle(pb, test_eq=False)
        self.assertMSONable(pb)

        self.assertTrue(pb.symbol == "Pb")
        self.assert_equal(pb.Z, 82.0)
        self.assert_equal(pb.Z_val, 14.0)
        self.assertTrue(pb.isnc)
        self.assertFalse(pb.ispaw)
        self.assert_equal(pb.l_max, 2)
        self.assert_equal(pb.l_local, 4)
        self.assertTrue(pb.supports_soc)


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
        #self.assertMSONable(table)

        selected = table.select_symbols("Si")
        assert len(selected) == len(table) and selected.__class__ is table.__class__

        with self.assertRaises(ValueError):
            table.pseudos_with_symbols("Si")
