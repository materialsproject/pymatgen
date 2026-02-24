from __future__ import annotations

import re

import pytest
from pytest import approx

from pymatgen.core import SETTINGS, Structure
from pymatgen.io.lobster.sets import LobsterSet
from pymatgen.io.vasp.sets import MODULE_DIR, BadInputSetWarning
from pymatgen.util.testing import FAKE_POTCAR_DIR, TEST_FILES_DIR, VASP_IN_DIR, MatSciTest

TEST_DIR = f"{TEST_FILES_DIR}/io/vasp"

pytest.MonkeyPatch().setitem(SETTINGS, "PMG_VASP_PSP_DIR", str(FAKE_POTCAR_DIR))

NO_PSP_DIR = SETTINGS.get("PMG_VASP_PSP_DIR") is None
skip_if_no_psp_dir = pytest.mark.skipif(NO_PSP_DIR, reason="PMG_VASP_PSP_DIR is not set")


class TestLobsterSet(MatSciTest):
    def setup_method(self):
        self.set = LobsterSet
        file_path = f"{VASP_IN_DIR}/POSCAR"
        file_path2 = f"{VASP_IN_DIR}/POSCAR.lobster.spin_DOS"
        self.struct = Structure.from_file(file_path)
        self.struct2 = Structure.from_file(file_path2)

        # test for different parameters!
        self.lobsterset1 = self.set(self.struct, isym=-1, ismear=-5)
        self.lobsterset2 = self.set(self.struct, isym=0, ismear=0)
        # only allow isym=-1 and isym=0
        with pytest.raises(
            ValueError,
            match=re.escape("Lobster cannot digest WAVEFUNCTIONS with symmetry. isym must be -1 or 0"),
        ):
            self.lobsterset_new = self.set(self.struct, isym=2, ismear=0)
        with pytest.raises(ValueError, match="Lobster usually works with ismear=-5 or ismear=0"):
            self.lobsterset_new = self.set(self.struct, isym=-1, ismear=2)
        # test if one can still hand over grid density of kpoints
        self.lobsterset3 = self.set(self.struct, isym=0, ismear=0, user_kpoints_settings={"grid_density": 6000})
        # check if users can overwrite settings in this class with the help of user_incar_settings
        self.lobsterset4 = self.set(self.struct, user_incar_settings={"ALGO": "Fast"})
        # use basis functions supplied by user
        self.lobsterset5 = self.set(
            self.struct,
            user_supplied_basis={"Fe": "3d 3p 4s", "P": "3p 3s", "O": "2p 2s"},
        )
        with pytest.raises(ValueError, match="There are no basis functions for the atom type O"):
            self.lobsterset6 = self.set(self.struct, user_supplied_basis={"Fe": "3d 3p 4s", "P": "3p 3s"}).incar
        self.lobsterset7 = self.set(
            self.struct,
            address_basis_file=f"{MODULE_DIR}/../lobster/lobster_basis/BASIS_PBE_54_standard.yaml",
        )
        with pytest.warns(BadInputSetWarning, match="Overriding the POTCAR"):
            self.lobsterset6 = self.set(self.struct)

        # test W_sw
        self.lobsterset8 = self.set(Structure.from_file(f"{TEST_FILES_DIR}/electronic_structure/cohp/POSCAR.W"))

        # test if potcar selection is consistent with PBE_54
        self.lobsterset9 = self.set(self.struct2)

    def test_incar(self):
        incar1 = self.lobsterset1.incar
        assert "NBANDS" in incar1
        assert incar1["NBANDS"] == 116
        assert incar1["NSW"] == 0
        assert incar1["ISMEAR"] == -5
        assert incar1["ISYM"] == -1
        assert incar1["ALGO"] == "Normal"
        assert incar1["EDIFF"] == approx(1e-6)
        incar2 = self.lobsterset2.incar
        assert incar2["ISYM"] == 0
        assert incar2["ISMEAR"] == 0
        incar4 = self.lobsterset4.incar
        assert incar4["ALGO"] == "Fast"

    def test_kpoints(self):
        kpoints1 = self.lobsterset1.kpoints
        assert kpoints1.comment.split()[5] == "6138"
        kpoints2 = self.lobsterset2.kpoints
        assert kpoints2.comment.split()[5] == "6138"
        kpoints3 = self.lobsterset3.kpoints
        assert kpoints3.comment.split()[5] == "6000"

    @skip_if_no_psp_dir
    def test_potcar(self):
        # PBE_54 is preferred at the moment
        functional, symbol = "PBE_54", "K_sv"
        assert self.lobsterset1.user_potcar_functional == functional
        # test if potcars selected are consistent with PBE_54
        assert self.lobsterset2.potcar.symbols == ["Fe_pv", "P", "O"]
        # test if error raised contains correct potcar symbol for K element as PBE_54 set
        with pytest.raises(
            FileNotFoundError,
            match=f"You do not have the right POTCAR with {functional=} and {symbol=}",
        ):
            _ = self.lobsterset9.potcar.symbols

    def test_as_from_dict(self):
        dict_here = self.lobsterset1.as_dict()

        lobsterset_new = self.set.from_dict(dict_here)
        # test relevant parts again
        incar1 = lobsterset_new.incar
        assert "NBANDS" in incar1
        assert incar1["NBANDS"] == 116
        assert incar1["NSW"] == 0
        assert incar1["NSW"] == 0
        assert incar1["ISMEAR"] == -5
        assert incar1["ISYM"] == -1
        assert incar1["ALGO"] == "Normal"
        kpoints1 = lobsterset_new.kpoints
        assert kpoints1.comment.split()[5] == "6138"
        assert lobsterset_new.user_potcar_functional == "PBE_54"
