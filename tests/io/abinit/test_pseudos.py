from __future__ import annotations

import os.path
import tarfile
from collections import defaultdict

import pytest
from monty.tempfile import ScratchDir
from pymatgen.io.abinit.pseudos import Pseudo, PseudoTable
from pymatgen.util.testing import TEST_FILES_DIR, PymatgenTest
from pytest import approx

TEST_DIR = f"{TEST_FILES_DIR}/io/abinit"


class TestPseudo(PymatgenTest):
    def setUp(self):
        nc_pseudo_fnames = defaultdict(list)
        nc_pseudo_fnames["Si"] = [f"{TEST_DIR}/{file}" for file in ("14si.pspnc", "14si.4.hgh", "14-Si.LDA.fhi")]

        self.nc_pseudos = defaultdict(list)
        self.paw_pseudos = defaultdict(list)

        for symbol, file_names in nc_pseudo_fnames.items():
            for file_name in file_names:
                _root, ext = os.path.splitext(file_name)
                pseudo = Pseudo.from_file(file_name)
                self.nc_pseudos[symbol].append(pseudo)

                # Save the pseudo as instance attribute whose name
                # is constructed with the rule: symbol_ppformat
                attribute = f"{symbol}_{ext[1:]}"
                if hasattr(self, attribute):
                    raise RuntimeError(f"self has already {attribute=}")

                setattr(self, attribute, pseudo)

    def test_nc_pseudos(self):
        """Test norm-conserving pseudopotentials."""
        for symbol, pseudos in self.nc_pseudos.items():
            for pseudo in pseudos:
                assert repr(pseudo)
                assert str(pseudo)
                assert pseudo.isnc
                assert not pseudo.ispaw
                assert pseudo.Z == 14
                assert pseudo.symbol == symbol
                assert pseudo.Z_val == 4
                assert pseudo.nlcc_radius >= 0.0

                # Test pickle
                self.serialize_with_pickle(pseudo)

                # Test MSONable
                self.assert_msonable(pseudo)

        # HGH pseudos
        pseudo = self.Si_hgh
        assert not pseudo.has_nlcc
        assert pseudo.l_max == 1
        assert pseudo.l_local == 0
        assert not pseudo.supports_soc
        assert self.Si_hgh.md5 is not None
        assert self.Si_hgh == self.Si_hgh

        # TM pseudos
        pseudo = self.Si_pspnc
        assert pseudo.has_nlcc
        assert pseudo.l_max == 2
        assert pseudo.l_local == 2
        assert not pseudo.supports_soc
        assert self.Si_hgh != self.Si_pspnc

        # FHI pseudos
        pseudo = self.Si_fhi
        assert not pseudo.has_nlcc
        assert pseudo.l_max == 3
        assert pseudo.l_local == 2
        assert not pseudo.supports_soc

        # Test PseudoTable
        table = PseudoTable(self.nc_pseudos["Si"])
        assert repr(table)
        assert str(table)
        assert table.allnc
        assert not table.allpaw
        assert table.is_complete
        assert len(table) == 3
        assert len(table[14]) == 3
        assert len(table.select_symbols("Si")) == 3
        assert table.zlist == [14]

        # Test pickle
        self.serialize_with_pickle(table, test_eq=False)

    def test_paw_pseudos(self):
        """Test 28ni.paw."""
        file_name = f"{TEST_DIR}/28ni.paw.tar.xz"
        symbol = "Ni"
        with ScratchDir(".") as tmp_dir, tarfile.open(file_name, mode="r:xz") as t:
            t.extractall(tmp_dir)
            path = os.path.join(tmp_dir, "28ni.paw")
            pseudo = Pseudo.from_file(path)

            assert repr(pseudo)
            assert str(pseudo)
            assert not pseudo.isnc
            assert pseudo.ispaw
            assert pseudo.Z == 28
            assert pseudo.symbol == symbol
            assert pseudo.Z_val == 18
            assert pseudo.paw_radius >= 0.0

            assert pseudo.l_max == 2
            assert pseudo.l_local == 0
            assert pseudo.supports_soc
            assert pseudo.md5 is not None

            # Test pickle
            self.serialize_with_pickle(pseudo)

            # Test MSONable
            self.assert_msonable(pseudo)

    def test_pawxml_pseudos(self):
        """Test O.GGA_PBE-JTH-paw.xml."""
        oxygen = Pseudo.from_file(f"{TEST_DIR}/O.GGA_PBE-JTH-paw.xml")
        assert repr(oxygen)
        assert str(oxygen)
        assert isinstance(oxygen.as_dict(), dict)

        assert oxygen.ispaw
        assert oxygen.symbol == "O"
        assert (oxygen.Z, oxygen.core, oxygen.valence) == (8, 2, 6), oxygen.Z_val == 6

        assert oxygen.xc.type == "GGA"
        assert oxygen.xc.name == "PBE"
        assert oxygen.supports_soc
        assert oxygen.md5 is not None
        assert oxygen.paw_radius == approx(1.4146523028)

        # Test pickle
        new_objs = self.serialize_with_pickle(oxygen)
        # Test MSONable
        self.assert_msonable(oxygen)

        for obj in new_objs:
            assert obj.ispaw
            assert obj.symbol == "O"
            assert (obj.Z, obj.core, obj.valence) == (8, 2, 6), obj.Z_val == 6

            assert obj.paw_radius == approx(1.4146523028)

    def test_oncvpsp_pseudo_sr(self):
        """Test the ONCVPSP Ge pseudo (scalar relativistic version)."""
        ger = Pseudo.from_file(f"{TEST_DIR}/ge.oncvpsp")
        assert repr(ger)
        assert str(ger)
        assert isinstance(ger.as_dict(), dict)
        ger.as_tmpfile()

        assert ger.symbol == "Ge"
        assert ger.Z == 32.0
        assert ger.Z_val == 4.0
        assert ger.isnc
        assert not ger.ispaw
        assert ger.l_max == 2
        assert ger.l_local == 4
        assert ger.rcore is None
        assert not ger.supports_soc

        # Data persistence
        self.serialize_with_pickle(ger)
        self.assert_msonable(ger)

    def test_oncvpsp_pseudo_fr(self):
        """Test the ONCVPSP Pb pseudo (relativistic version with SO)."""
        pb = Pseudo.from_file(f"{TEST_DIR}/Pb-d-3_r.psp8")
        repr(pb)
        str(pb)

        # Data persistence
        self.serialize_with_pickle(pb)
        self.assert_msonable(pb)

        assert pb.symbol == "Pb"
        assert pb.Z == 82.0
        assert pb.Z_val == 14.0
        assert pb.isnc
        assert not pb.ispaw
        assert pb.l_max == 2
        assert pb.l_local == 4
        assert pb.supports_soc


class TestPseudoTable(PymatgenTest):
    def test_methods(self):
        """Test PseudoTable methods."""
        table = PseudoTable([f"{TEST_DIR}/{file}" for file in ("14si.pspnc", "14si.4.hgh", "14-Si.LDA.fhi")])
        assert str(table)
        assert len(table) == 3
        for pseudo in table:
            assert pseudo.isnc
        assert table.allnc
        assert not table.allpaw
        assert table.zlist == [14]

        # Data persistence
        self.serialize_with_pickle(table, test_eq=False)

        dct = table.as_dict()
        PseudoTable.from_dict(dct)
        self.assert_msonable(table)

        selected = table.select_symbols("Si")
        assert len(selected) == len(table)
        assert selected.__class__ is table.__class__

        with pytest.raises(ValueError, match=r"Found multiple occurrences of symbol\(s\) Si"):
            table.pseudos_with_symbols("Si")
