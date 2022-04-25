# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.


import collections
import os.path

from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable
from pymatgen.util.testing import PymatgenTest

_test_dir = os.path.join(PymatgenTest.TEST_FILES_DIR, "abinit")


def ref_file(filename):
    return os.path.join(_test_dir, filename)


def ref_files(*filenames):
    return list(map(ref_file, filenames))


class PseudoTestCase(PymatgenTest):
    def setUp(self):
        nc_pseudo_fnames = collections.defaultdict(list)
        nc_pseudo_fnames["Si"] = ref_files("14si.pspnc", "14si.4.hgh", "14-Si.LDA.fhi")

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
                    raise RuntimeError(f"self has already the attribute {attr_name}")

                setattr(self, attr_name, pseudo)

    def test_nc_pseudos(self):
        """Test norm-conserving pseudopotentials"""
        for symbol, pseudos in self.nc_pseudos.items():
            for pseudo in pseudos:
                assert repr(pseudo)
                assert str(pseudo)
                self.assertTrue(pseudo.isnc)
                self.assertFalse(pseudo.ispaw)
                self.assertEqual(pseudo.Z, 14)
                self.assertEqual(pseudo.symbol, symbol)
                self.assertEqual(pseudo.Z_val, 4)
                self.assertGreaterEqual(pseudo.nlcc_radius, 0.0)

                # Test pickle
                self.serialize_with_pickle(pseudo, test_eq=False)

                # Test MSONable
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
        assert repr(table)
        assert str(table)
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
        assert repr(oxygen)
        assert str(oxygen)
        assert isinstance(oxygen.as_dict(), dict)

        self.assertTrue(oxygen.ispaw)
        self.assertTrue(
            oxygen.symbol == "O" and (oxygen.Z, oxygen.core, oxygen.valence) == (8, 2, 6),
            oxygen.Z_val == 6,
        )

        assert oxygen.xc.type == "GGA" and oxygen.xc.name == "PBE"
        assert oxygen.supports_soc
        assert oxygen.md5 is not None
        self.assertAlmostEqual(oxygen.paw_radius, 1.4146523028)

        # Test pickle
        new_objs = self.serialize_with_pickle(oxygen, test_eq=False)
        # Test MSONable
        self.assertMSONable(oxygen)

        for o in new_objs:
            self.assertTrue(o.ispaw)
            self.assertTrue(
                o.symbol == "O" and (o.Z, o.core, o.valence) == (8, 2, 6),
                o.Z_val == 6,
            )

            self.assertAlmostEqual(o.paw_radius, 1.4146523028)

    def test_oncvpsp_pseudo_sr(self):
        """
        Test the ONCVPSP Ge pseudo (scalar relativistic version).
        """
        ger = Pseudo.from_file(ref_file("ge.oncvpsp"))
        assert repr(ger)
        assert str(ger)
        assert isinstance(ger.as_dict(), dict)
        ger.as_tmpfile()

        self.assertTrue(ger.symbol == "Ge")
        self.assertEqual(ger.Z, 32.0)
        self.assertEqual(ger.Z_val, 4.0)
        self.assertTrue(ger.isnc)
        self.assertFalse(ger.ispaw)
        self.assertEqual(ger.l_max, 2)
        self.assertEqual(ger.l_local, 4)
        self.assertEqual(ger.rcore, None)
        assert not ger.supports_soc

        # Data persistence
        self.serialize_with_pickle(ger, test_eq=False)
        self.assertMSONable(ger)

    def test_oncvpsp_pseudo_fr(self):
        """
        Test the ONCVPSP Pb pseudo (relativistic version with SO).
        """
        pb = Pseudo.from_file(ref_file("Pb-d-3_r.psp8"))
        repr(pb)
        str(pb)

        # Data persistence
        self.serialize_with_pickle(pb, test_eq=False)
        self.assertMSONable(pb)

        self.assertTrue(pb.symbol == "Pb")
        self.assertEqual(pb.Z, 82.0)
        self.assertEqual(pb.Z_val, 14.0)
        self.assertTrue(pb.isnc)
        self.assertFalse(pb.ispaw)
        self.assertEqual(pb.l_max, 2)
        self.assertEqual(pb.l_local, 4)
        self.assertTrue(pb.supports_soc)


class PseudoTableTest(PymatgenTest):
    def test_methods(self):
        """Test PseudoTable methods"""
        table = PseudoTable(ref_files("14si.pspnc", "14si.4.hgh", "14-Si.LDA.fhi"))
        assert str(table)
        assert len(table) == 3
        for pseudo in table:
            assert pseudo.isnc
        assert table.allnc and not table.allpaw
        assert table.zlist == [14]

        # Data persistence
        self.serialize_with_pickle(table, test_eq=False)

        d = table.as_dict()
        PseudoTable.from_dict(d)
        self.assertMSONable(table)

        selected = table.select_symbols("Si")
        assert len(selected) == len(table) and selected.__class__ is table.__class__

        with self.assertRaises(ValueError):
            table.pseudos_with_symbols("Si")
