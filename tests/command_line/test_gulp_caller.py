"""
Created on Jan 22, 2013.

@author: Bharat Medasani
"""

from __future__ import annotations

import os
import sys
import unittest
from shutil import which

import pytest

from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.command_line.gulp_caller import (
    BuckinghamPotential,
    GulpCaller,
    GulpError,
    GulpIO,
    get_energy_buckingham,
    get_energy_relax_structure_buckingham,
    get_energy_tersoff,
)
from pymatgen.core.structure import Structure
from pymatgen.util.testing import TEST_FILES_DIR

gulp_present = which("gulp") and os.getenv("GULP_LIB") and ("win" not in sys.platform)
# disable gulp tests for now. Right now, it is compiled against libgfortran3, which is no longer supported in the new
# Ubuntu 20.04.
gulp_present = False


@unittest.skipIf(not gulp_present, "gulp not present.")
class TestGulpCaller(unittest.TestCase):
    def test_run(self):
        mgo_lattice = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg"] * 4 + ["O"] * 4
        mgo_frac_cord = [
            [0, 0, 0],
            [0.5, 0.5, 0],
            [0.5, 0, 0.5],
            [0, 0.5, 0.5],
            [0.5, 0, 0],
            [0, 0.5, 0],
            [0, 0, 0.5],
            [0.5, 0.5, 0.5],
        ]
        mgo_uc = Structure(mgo_lattice, mgo_specie, mgo_frac_cord, validate_proximity=True, to_unit_cell=True)
        gio = GulpIO()
        gin = gio.keyword_line("optimise", "conp")
        gin += gio.structure_lines(mgo_uc, symm_flg=False)
        # gin += self.gc.gulp_lib('catlow.lib')
        gin += "species\nMg    core  2.00000\nO core  0.86902\nO shel -2.86902\n"
        gin += "buck\n"
        gin += "Mg core O shel   946.627 0.31813  0.00000 0.0 10.0\n"
        gin += "O  shel O shel 22764.000 0.14900 27.87900 0.0 12.0\n"
        gc = GulpCaller()

        """Some inherent checks are in the run_gulp function itself.
        They should be sufficient for raising errors."""
        gc.run(gin)

    def test_decimal(self):
        struct = Structure.from_str(
            """Mg2 Al4 O8
        1.0
        5.003532 0.000000 2.888790
        1.667844 4.717375 2.888790
        0.000000 0.000000 5.777581
        O Mg Al
        8 2 4
        direct
        0.736371 0.736371 0.736371 O
        0.263629 0.263629 0.709114 O
        0.263629 0.709114 0.263629 O
        0.709114 0.263629 0.263629 O
        0.736371 0.290886 0.736371 O
        0.290886 0.736371 0.736371 O
        0.263629 0.263629 0.263629 O
        0.736371 0.736371 0.290886 O
        0.125000 0.125000 0.125000 Mg
        0.875000 0.875000 0.875000 Mg
        0.500000 0.500000 0.000000 Al
        0.500000 0.500000 0.500000 Al
        0.000000 0.500000 0.500000 Al
        0.500000 0.000000 0.500000 Al""",
            fmt="poscar",
        )

        _ = BuckinghamPotential(bush_lewis_flag="bush")
        gio = GulpIO()
        buckingham_input = gio.buckingham_input(struct, ["relax conp"])
        caller = GulpCaller()
        caller.run(buckingham_input)


@unittest.skipIf(not gulp_present, "gulp not present.")
class TestGulpIO(unittest.TestCase):
    def setUp(self):
        self.structure = Structure.from_file(f"{TEST_FILES_DIR}/POSCAR.Al12O18")
        self.gio = GulpIO()

    def test_keyword_line_with_correct_keywords(self):
        kw = ("defect", "property")
        inp_str = self.gio.keyword_line(*kw)
        for word in kw:
            assert word in inp_str

    def test_structure_lines_default_options(self):
        inp_str = self.gio.structure_lines(self.structure)
        assert "cell" in inp_str
        assert "frac" in inp_str
        assert "space" in inp_str

    def test_structure_lines_no_unitcell(self):
        inp_str = self.gio.structure_lines(self.structure, cell_flg=False)
        assert "cell" not in inp_str

    def test_structure_lines_no_frac_coords(self):
        inp_str = self.gio.structure_lines(self.structure, cell_flg=False, frac_flg=False)
        assert "cell" not in inp_str
        assert "cart" in inp_str

    @unittest.skip("Not Implemented yet")
    def test_specie_potential(self):
        pass

    @unittest.expectedFailure
    def test_library_line_explicit_path(self):
        gin = self.gio.library_line("/Users/mbkumar/Research/Defects/GulpExe/Libraries/catlow.lib")
        assert "lib" in gin

    def test_library_line_wrong_file(self):
        with pytest.raises(GulpError, match="GULP library not found"):
            self.gio.library_line("temp_to_fail.lib")

    def test_buckingham_potential(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg", "O"] * 4
        mgo_frac_cord = [
            [0, 0, 0],
            [0.5, 0, 0],
            [0.5, 0.5, 0],
            [0, 0.5, 0],
            [0.5, 0, 0.5],
            [0, 0, 0.5],
            [0, 0.5, 0.5],
            [0.5, 0.5, 0.5],
        ]
        mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, validate_proximity=True, to_unit_cell=True)
        gin = self.gio.buckingham_potential(mgo_uc)
        assert "specie" in gin
        assert "buck" in gin
        assert "spring" in gin
        assert "Mg core" in gin
        assert "O  core" in gin
        assert "O  shel" in gin

        gin = self.gio.buckingham_potential(self.structure)
        assert "specie" in gin
        assert "buck" in gin
        assert "spring" in gin

    def test_buckingham_input(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg", "O"] * 4
        mgo_frac_cord = [
            [0, 0, 0],
            [0.5, 0, 0],
            [0.5, 0.5, 0],
            [0, 0.5, 0],
            [0.5, 0, 0.5],
            [0, 0, 0.5],
            [0, 0.5, 0.5],
            [0.5, 0.5, 0.5],
        ]
        mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, validate_proximity=True, to_unit_cell=True)
        gin = self.gio.buckingham_input(mgo_uc, keywords=("optimise", "conp"))
        assert "optimise" in gin
        assert "cell" in gin
        assert "specie" in gin
        assert "buck" in gin
        assert "spring" in gin
        assert "Mg core" in gin
        assert "O  core" in gin
        assert "O  shel" in gin

    # Improve the test
    def test_tersoff_potential(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg", "O"] * 4
        mgo_frac_cord = [
            [0, 0, 0],
            [0.5, 0, 0],
            [0.5, 0.5, 0],
            [0, 0.5, 0],
            [0.5, 0, 0.5],
            [0, 0, 0.5],
            [0, 0.5, 0.5],
            [0.5, 0.5, 0.5],
        ]
        mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, validate_proximity=True, to_unit_cell=True)
        gin = self.gio.tersoff_potential(mgo_uc)
        assert "specie" in gin
        assert "Mg core" in gin

    def test_get_energy(self):
        # Output string obtained from running GULP on a terminal
        out_str = """  Components of energy :
--------------------------------------------------------------------------------
  Interatomic potentials     =           5.61135426 eV
  Monopole - monopole (real) =          -4.34238722 eV
  Monopole - monopole (recip)=         -43.45344934 eV
  Monopole - monopole (total)=         -47.79583656 eV
--------------------------------------------------------------------------------
  Total lattice energy :
    Primitive unit cell      =         -42.18448230 eV
    Non-primitive unit cell  =        -168.73792920 eV
--------------------------------------------------------------------------------
  Total lattice energy (in kJmol-1):
    Primitive unit cell      =           -4070.1577 kJ/(mole unit cells)
    Non-primitive unit cell  =          -16280.6308 kJ/(mole unit cells)
--------------------------------------------------------------------------------
  Components of energy :

--------------------------------------------------------------------------------
  Interatomic potentials     =           6.79846039 eV
  Monopole - monopole (real) =          -4.45761741 eV
  Monopole - monopole (recip)=         -44.60653603 eV
  Monopole - monopole (total)=         -49.06415344 eV
--------------------------------------------------------------------------------
  Total lattice energy :
    Primitive unit cell      =         -42.26569304 eV
    Non-primitive unit cell  =        -169.06277218 eV
--------------------------------------------------------------------------------
  Total lattice energy (in kJmol-1):
    Primitive unit cell      =           -4077.9933 kJ/(mole unit cells)
    Non-primitive unit cell  =          -16311.9732 kJ/(mole unit cells)
--------------------------------------------------------------------------------"""
        energy = self.gio.get_energy(out_str)
        assert energy == -169.06277218

    def test_get_relaxed_structure(self):
        # Output string obtained from running GULP on a terminal

        with open(f"{TEST_FILES_DIR}/example21.gout") as fp:
            out_str = fp.read()
        struct = self.gio.get_relaxed_structure(out_str)
        assert isinstance(struct, Structure)
        assert len(struct) == 8
        assert struct.lattice.a == 4.212
        assert struct.lattice.alpha == 90

    @unittest.skip("Test later")
    def test_tersoff_input(self):
        self.gio.tersoff_input(self.structure)


@unittest.skipIf(not gulp_present, "gulp not present.")
class TestGlobalFunctions(unittest.TestCase):
    def setUp(self):
        mgo_latt = [[4.212, 0, 0], [0, 4.212, 0], [0, 0, 4.212]]
        mgo_specie = ["Mg", "O"] * 4
        mgo_frac_cord = [
            [0, 0, 0],
            [0.5, 0, 0],
            [0.5, 0.5, 0],
            [0, 0.5, 0],
            [0.5, 0, 0.5],
            [0, 0, 0.5],
            [0, 0.5, 0.5],
            [0.5, 0.5, 0.5],
        ]
        self.mgo_uc = Structure(mgo_latt, mgo_specie, mgo_frac_cord, validate_proximity=True, to_unit_cell=True)
        bv = BVAnalyzer()
        val = bv.get_valences(self.mgo_uc)
        el = [site.species_string for site in self.mgo_uc]
        self.val_dict = dict(zip(el, val))

    def test_get_energy_tersoff(self):
        structure = Structure.from_file(f"{TEST_FILES_DIR}/POSCAR.Al12O18")
        energy = get_energy_tersoff(structure)
        assert isinstance(energy, float)

    def test_get_energy_buckingham(self):
        energy = get_energy_buckingham(self.mgo_uc)
        assert isinstance(energy, float)
        # test with vacancy structure
        del self.mgo_uc[0]
        energy = get_energy_buckingham(
            self.mgo_uc,
            keywords=("qok", "optimise", "conp"),
            valence_dict=self.val_dict,
        )
        assert isinstance(energy, float)

    def test_get_energy_relax_structure_buckingham(self):
        energy, struct = get_energy_relax_structure_buckingham(self.mgo_uc)
        assert isinstance(energy, float)
        assert isinstance(struct, Structure)
        site_len = len(struct)
        assert site_len == len(self.mgo_uc)


@unittest.skipIf(not gulp_present, "gulp not present.")
class TestBuckinghamPotentialLewis(unittest.TestCase):
    def setUp(self):
        self.bpl = BuckinghamPotential("lewis")

    def test_existing_element(self):
        assert "Sc_2+" in self.bpl.pot_dict
        assert "Sc_2+" in self.bpl.species_dict
        assert "O" in self.bpl.pot_dict
        assert "O_core" in self.bpl.species_dict
        assert "O_shel" in self.bpl.species_dict

    def test_non_existing_element(self):
        assert "Li_1+" not in self.bpl.pot_dict
        assert "Li_1+" not in self.bpl.species_dict

    def test_element_different_valence(self):
        assert "Sc_4+" not in self.bpl.species_dict

    def test_values(self):
        assert self.bpl.species_dict["Sc_2+"] != ""
        assert self.bpl.pot_dict["Sc_2+"] != ""

    def test_spring(self):
        assert "Li" not in self.bpl.spring_dict
        assert self.bpl.spring_dict["O"] != ""


@unittest.skipIf(not gulp_present, "gulp not present.")
class TestBuckinghamPotentialBush(unittest.TestCase):
    def setUp(self):
        self.bpb = BuckinghamPotential("bush")

    def test_existing_element(self):
        assert "Li" in self.bpb.pot_dict
        assert "Li" in self.bpb.species_dict
        assert "O" in self.bpb.pot_dict
        assert "O" in self.bpb.species_dict

    def test_non_existing_element(self):
        assert "Mn" not in self.bpb.pot_dict
        assert "Mn" not in self.bpb.species_dict

    def test_element_different_valence(self):
        assert self.bpb.species_dict["Li"]["oxi"] != 2

    def test_spring(self):
        assert self.bpb.spring_dict["Li"] == ""
        assert self.bpb.spring_dict["O"] != ""
