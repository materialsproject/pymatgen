from __future__ import annotations

import os
import pickle
from pathlib import Path

import numpy as np
import pytest
import scipy.constants as const
from monty.io import zopen
from monty.serialization import loadfn
from pytest import approx

import pymatgen
from pymatgen.core import SETTINGS
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.core import Magmom
from pymatgen.io.vasp.inputs import (
    BadIncarWarning,
    Incar,
    Kpoints,
    Poscar,
    Potcar,
    PotcarSingle,
    UnknownPotcarWarning,
    VaspInput,
)
from pymatgen.util.testing import TEST_FILES_DIR


class TestPoscar:
    def test_init(self):
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        poscar = Poscar.from_file(filepath, check_for_POTCAR=False)
        comp = poscar.structure.composition
        assert comp == Composition("Fe4P4O16")

        # VASP 4 type with symbols at the end.
        poscar_string = """Test1
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 F
"""
        poscar = Poscar.from_str(poscar_string)
        assert poscar.structure.composition == Composition("SiF")

        poscar_string = ""
        with pytest.raises(ValueError, match="Empty POSCAR"):
            Poscar.from_str(poscar_string)

        # VASP 4 style file with default names, i.e. no element symbol found.
        poscar_string = """Test2
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
direct
0.000000 0.000000 0.000000
0.750000 0.500000 0.750000
"""
        poscar = Poscar.from_str(poscar_string)
        assert poscar.structure.composition == Composition("HHe")
        # VASP 4 style file with default names, i.e. no element symbol found.
        poscar_string = """Test3
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
Selective dynamics
direct
0.000000 0.000000 0.000000 T T T Si
0.750000 0.500000 0.750000 F F F O
"""
        poscar = Poscar.from_str(poscar_string)
        selective_dynamics = [list(x) for x in poscar.selective_dynamics]

        assert selective_dynamics == [[True, True, True], [False, False, False]]
        self.selective_poscar = poscar

    def test_from_file(self):
        filepath = f"{TEST_FILES_DIR}/POSCAR.symbols_natoms_multilines"
        poscar = Poscar.from_file(filepath, check_for_POTCAR=False, read_velocities=False)
        ordered_expected_elements = [
            "Fe",
            "Cr",
            "Fe",
            "Fe",
            "Cr",
            "Cr",
            "Cr",
            "Cr",
            "Fe",
            "Fe",
            "Cr",
            "Fe",
            "Cr",
            "Fe",
            "Fe",
            "Cr",
            "Fe",
            "Cr",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Cr",
            "Fe",
            "Ni",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Cr",
            "Cr",
            "Cr",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Cr",
            "Fe",
            "Fe",
            "Ni",
            "Fe",
            "Fe",
            "Fe",
            "Cr",
            "Cr",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
            "Fe",
        ]
        assert [site.specie.symbol for site in poscar.structure] == ordered_expected_elements

    def test_to_from_dict(self):
        poscar_string = """Test3
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
Selective dynamics
direct
0.000000 0.000000 0.000000 T T T Si
0.750000 0.500000 0.750000 F F F O
"""
        poscar = Poscar.from_str(poscar_string)
        d = poscar.as_dict()
        poscar2 = Poscar.from_dict(d)
        assert poscar2.comment == "Test3"
        assert all(poscar2.selective_dynamics[0])
        assert not all(poscar2.selective_dynamics[1])

    def test_cart_scale(self):
        poscar_string = """Test1
1.1
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si F
1 1
cart
0.000000   0.00000000   0.00000000
3.840198   1.50000000   2.35163175
"""
        p = Poscar.from_str(poscar_string)
        site = p.structure[1]
        assert np.allclose(site.coords, np.array([3.840198, 1.5, 2.35163175]) * 1.1)

    def test_significant_figures(self):
        si = 14
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])

        # Silicon structure for testing.
        latt = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(latt, [si, si], coords)
        poscar = Poscar(struct)
        expected_str = """Si2
1.0
   3.84    0.00    0.00
   1.92    3.33    0.00
   0.00   -2.22    3.14
Si
2
direct
   0.00    0.00    0.00 Si
   0.75    0.50    0.75 Si
"""

        actual_str = poscar.get_string(significant_figures=2)
        assert actual_str == expected_str, "Wrong POSCAR output!"

    def test_str(self):
        si = 14
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])

        # Silicon structure for testing.
        latt = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(latt, [si, si], coords)
        poscar = Poscar(struct)
        expected_str = """Si2
1.0
   3.8401979336999998    0.0000000000000000    0.0000000000000000
   1.9200989667999999    3.3257101909000002    0.0000000000000000
   0.0000000000000000   -2.2171384942999999    3.1355090603000000
Si
2
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Si
   0.7500000000000000    0.5000000000000000    0.7500000000000000 Si
"""

        assert str(poscar) == expected_str, "Wrong POSCAR output!"

        # VASP 4 type with symbols at the end.
        poscar_string = """Test1
1.0
-3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 F
"""

        expected = """Test1
1.0
   3.8401980000000000   -0.0000000000000000   -0.0000000000000000
  -1.9200990000000000   -3.3257099999999999   -0.0000000000000000
  -0.0000000000000000    2.2171379999999998   -3.1355089999999999
Si F
1 1
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Si
   0.7500000000000000    0.5000000000000000    0.7500000000000000 F
"""
        poscar = Poscar.from_str(poscar_string)
        assert str(poscar) == expected

    def test_from_md_run(self):
        # Parsing from an MD type run with velocities and predictor corrector data
        p = Poscar.from_file(f"{TEST_FILES_DIR}/CONTCAR.MD", check_for_POTCAR=False)
        assert np.sum(np.array(p.velocities)) == approx(0.0065417961324)
        assert p.predictor_corrector[0][0][0] == 0.33387820e00
        assert p.predictor_corrector[0][1][1] == -0.10583589e-02

    def test_write_md_poscar(self):
        # Parsing from an MD type run with velocities and predictor corrector data
        # And writing a new POSCAR from the new structure
        p = Poscar.from_file(f"{TEST_FILES_DIR}/CONTCAR.MD", check_for_POTCAR=False)

        path = Path("POSCAR.testing.md")
        p.write_file(path)
        p3 = Poscar.from_file(path)

        assert np.allclose(p.structure.lattice.abc, p3.structure.lattice.abc, 5)
        assert np.allclose(p.velocities, p3.velocities, 5)
        assert np.allclose(p.predictor_corrector, p3.predictor_corrector, 5)
        assert p.predictor_corrector_preamble == p3.predictor_corrector_preamble
        path.unlink()

    def test_setattr(self):
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        poscar = Poscar.from_file(filepath, check_for_POTCAR=False)
        with pytest.raises(ValueError, match="velocities array must be same length as the structure"):
            poscar.velocities = [[0, 0, 0]]
        poscar.selective_dynamics = np.array([[True, False, False]] * 24)
        expected = """
        Fe4P4O16
1.0
  10.4117668699494264    0.0000000000000000    0.0000000000000000
   0.0000000000000000    6.0671718799705294    0.0000000000000000
   0.0000000000000000    0.0000000000000000    4.7594895399768813
Fe P O
4 4 16
Selective dynamics
direct
   0.2187282200000000    0.7500000000000000    0.4748671100000000 T F F Fe
   0.2812717800000000    0.2500000000000000    0.9748671100000000 T F F Fe
   0.7187282200000000    0.7500000000000000    0.0251328900000000 T F F Fe
   0.7812717800000000    0.2500000000000000    0.5251328900000000 T F F Fe
   0.0946130900000000    0.2500000000000000    0.4182432700000000 T F F P
   0.4053869100000000    0.7500000000000000    0.9182432699999999 T F F P
   0.5946130900000000    0.2500000000000000    0.0817567300000000 T F F P
   0.9053869100000000    0.7500000000000000    0.5817567300000001 T F F P
   0.0433723100000000    0.7500000000000000    0.7071376700000001 T F F O
   0.0966424400000000    0.2500000000000000    0.7413203500000000 T F F O
   0.1657097400000000    0.0460723300000000    0.2853839400000000 T F F O
   0.1657097400000000    0.4539276700000000    0.2853839400000000 T F F O
   0.3342902600000000    0.5460723300000000    0.7853839400000000 T F F O
   0.3342902600000000    0.9539276700000000    0.7853839400000000 T F F O
   0.4033575600000000    0.7500000000000000    0.2413203500000000 T F F O
   0.4566276900000000    0.2500000000000000    0.2071376700000000 T F F O
   0.5433723100000000    0.7500000000000000    0.7928623299999999 T F F O
   0.5966424400000000    0.2500000000000000    0.7586796500000000 T F F O
   0.6657097400000000    0.0460723300000000    0.2146160600000000 T F F O
   0.6657097400000000    0.4539276700000000    0.2146160600000000 T F F O
   0.8342902600000000    0.5460723300000000    0.7146160600000000 T F F O
   0.8342902600000000    0.9539276700000000    0.7146160600000000 T F F O
   0.9033575600000000    0.7500000000000000    0.2586796500000000 T F F O
   0.9566276900000000    0.2500000000000000    0.2928623300000000 T F F O"""
        assert str(poscar).strip() == expected.strip()
        poscar.velocities = np.ones((24, 3))
        assert "velocities" in poscar.structure.site_properties

    def test_velocities(self):
        si = 14
        coords = []
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])

        # Silicon structure for testing.
        latt = [
            [3.8401979337, 0.00, 0.00],
            [1.9200989668, 3.3257101909, 0.00],
            [0.00, -2.2171384943, 3.1355090603],
        ]
        struct = Structure(latt, [si, si], coords)
        poscar = Poscar(struct)
        poscar.set_temperature(900)

        v = np.array(poscar.velocities)

        for x in np.sum(v, axis=0):
            assert x == approx(0, abs=1e-7)

        temperature = struct[0].specie.atomic_mass.to("kg") * np.sum(v**2) / (3 * const.k) * 1e10
        assert temperature == approx(900, abs=1e-4), "Temperature instantiated incorrectly"

        poscar.set_temperature(700)
        v = np.array(poscar.velocities)
        for x in np.sum(v, axis=0):
            assert x == approx(0, abs=1e-7), "Velocities initialized with a net momentum"

        temperature = struct[0].specie.atomic_mass.to("kg") * np.sum(v**2) / (3 * const.k) * 1e10
        assert temperature == approx(700, abs=1e-4), "Temperature instantiated incorrectly"

    def test_write(self):
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        poscar = Poscar.from_file(filepath)
        tempfname = Path("POSCAR.testing")
        poscar.write_file(tempfname)
        p = Poscar.from_file(tempfname)
        assert np.allclose(poscar.structure.lattice.abc, p.structure.lattice.abc, 5)
        tempfname.unlink()

    def test_selective_dynamics(self):
        filepath = f"{TEST_FILES_DIR}/POSCAR.Fe3O4"
        poscar = Poscar.from_file(filepath)
        structure = poscar.structure

        # Fix bottom half
        fixed_indices = structure.frac_coords[:, 2] >= 0.5

        poscar = Poscar(structure, selective_dynamics=np.tile(fixed_indices.reshape(-1, 1), [1, 3]))
        selective_dynamics = [list(x) for x in poscar.selective_dynamics]

        assert selective_dynamics == [
            [True, True, True],
            [False, False, False],
            [False, False, False],
            [True, True, True],
            [False, False, False],
            [True, True, True],
            [True, True, True],
            [False, False, False],
            [False, False, False],
            [True, True, True],
            [True, True, True],
            [False, False, False],
            [True, True, True],
            [False, False, False],
        ]


class TestIncar:
    def setup(self):
        file_name = f"{TEST_FILES_DIR}/INCAR"
        self.incar = Incar.from_file(file_name)

    def test_init(self):
        incar = self.incar
        incar["LDAU"] = "T"
        assert incar["ALGO"] == "Damped", "Wrong Algo"
        assert float(incar["EDIFF"]) == 1e-4, "Wrong EDIFF"
        assert isinstance(incar["LORBIT"], int)

    def test_diff(self):
        filepath1 = f"{TEST_FILES_DIR}/INCAR"
        incar1 = Incar.from_file(filepath1)
        filepath2 = f"{TEST_FILES_DIR}/INCAR.2"
        incar2 = Incar.from_file(filepath2)
        incar3 = Incar.from_file(filepath2)
        assert incar1.diff(incar2) == {
            "Different": {
                "NELM": {"INCAR1": None, "INCAR2": 100},
                "ISPIND": {"INCAR1": 2, "INCAR2": None},
                "LWAVE": {"INCAR1": True, "INCAR2": False},
                "LDAUPRINT": {"INCAR1": None, "INCAR2": 1},
                "MAGMOM": {
                    "INCAR1": [
                        6,
                        -6,
                        -6,
                        6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                    ],
                    "INCAR2": None,
                },
                "NELMIN": {"INCAR1": None, "INCAR2": 3},
                "ENCUTFOCK": {"INCAR1": 0.0, "INCAR2": None},
                "HFSCREEN": {"INCAR1": 0.207, "INCAR2": None},
                "LSCALU": {"INCAR1": False, "INCAR2": None},
                "ENCUT": {"INCAR1": 500, "INCAR2": None},
                "NSIM": {"INCAR1": 1, "INCAR2": None},
                "ICHARG": {"INCAR1": None, "INCAR2": 1},
                "NSW": {"INCAR1": 99, "INCAR2": 51},
                "NKRED": {"INCAR1": 2, "INCAR2": None},
                "NUPDOWN": {"INCAR1": 0, "INCAR2": None},
                "LCHARG": {"INCAR1": True, "INCAR2": None},
                "LPLANE": {"INCAR1": True, "INCAR2": None},
                "ISMEAR": {"INCAR1": 0, "INCAR2": -5},
                "NPAR": {"INCAR1": 8, "INCAR2": 1},
                "SYSTEM": {
                    "INCAR1": "Id=[0] dblock_code=[97763-icsd] formula=[li mn (p o4)] sg_name=[p n m a]",
                    "INCAR2": "Id=[91090] dblock_code=[20070929235612linio-59.53134651-vasp] formula=[li3 ni3 o6] "
                    "sg_name=[r-3m]",
                },
                "ALGO": {"INCAR1": "Damped", "INCAR2": "Fast"},
                "LHFCALC": {"INCAR1": True, "INCAR2": None},
                "TIME": {"INCAR1": 0.4, "INCAR2": None},
            },
            "Same": {
                "IBRION": 2,
                "PREC": "Accurate",
                "ISIF": 3,
                "LMAXMIX": 4,
                "LREAL": "Auto",
                "ISPIN": 2,
                "EDIFF": 0.0001,
                "LORBIT": 11,
                "SIGMA": 0.05,
            },
        }

        assert incar1.diff(incar3) == {
            "Different": {
                "NELM": {"INCAR1": None, "INCAR2": 100},
                "ISPIND": {"INCAR1": 2, "INCAR2": None},
                "LWAVE": {"INCAR1": True, "INCAR2": False},
                "LDAUPRINT": {"INCAR1": None, "INCAR2": 1},
                "MAGMOM": {
                    "INCAR1": [
                        6,
                        -6,
                        -6,
                        6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                        0.6,
                    ],
                    "INCAR2": None,
                },
                "NELMIN": {"INCAR1": None, "INCAR2": 3},
                "ENCUTFOCK": {"INCAR1": 0.0, "INCAR2": None},
                "HFSCREEN": {"INCAR1": 0.207, "INCAR2": None},
                "LSCALU": {"INCAR1": False, "INCAR2": None},
                "ENCUT": {"INCAR1": 500, "INCAR2": None},
                "NSIM": {"INCAR1": 1, "INCAR2": None},
                "ICHARG": {"INCAR1": None, "INCAR2": 1},
                "NSW": {"INCAR1": 99, "INCAR2": 51},
                "NKRED": {"INCAR1": 2, "INCAR2": None},
                "NUPDOWN": {"INCAR1": 0, "INCAR2": None},
                "LCHARG": {"INCAR1": True, "INCAR2": None},
                "LPLANE": {"INCAR1": True, "INCAR2": None},
                "ISMEAR": {"INCAR1": 0, "INCAR2": -5},
                "NPAR": {"INCAR1": 8, "INCAR2": 1},
                "SYSTEM": {
                    "INCAR1": "Id=[0] dblock_code=[97763-icsd] formula=[li mn (p o4)] sg_name=[p n m a]",
                    "INCAR2": "Id=[91090] dblock_code=[20070929235612linio-59.53134651-vasp] formula=[li3 ni3 o6] "
                    "sg_name=[r-3m]",
                },
                "ALGO": {"INCAR1": "Damped", "INCAR2": "Fast"},
                "LHFCALC": {"INCAR1": True, "INCAR2": None},
                "TIME": {"INCAR1": 0.4, "INCAR2": None},
            },
            "Same": {
                "IBRION": 2,
                "PREC": "Accurate",
                "ISIF": 3,
                "LMAXMIX": 4,
                "LREAL": "Auto",
                "ISPIN": 2,
                "EDIFF": 0.0001,
                "LORBIT": 11,
                "SIGMA": 0.05,
            },
        }

    def test_as_dict_and_from_dict(self):
        d = self.incar.as_dict()
        incar2 = Incar.from_dict(d)
        assert self.incar == incar2
        d["MAGMOM"] = [Magmom([1, 2, 3]).as_dict()]
        incar3 = Incar.from_dict(d)
        assert incar3["MAGMOM"] == [Magmom([1, 2, 3])]

    def test_write(self):
        tempfname = Path("INCAR.testing")
        self.incar.write_file(tempfname)
        i = Incar.from_file(tempfname)
        assert i == self.incar
        tempfname.unlink()

    def test_get_string(self):
        s = self.incar.get_string(pretty=True, sort_keys=True)
        expected = """ALGO       =  Damped
EDIFF      =  0.0001
ENCUT      =  500
ENCUTFOCK  =  0.0
HFSCREEN   =  0.207
IBRION     =  2
ISIF       =  3
ISMEAR     =  0
ISPIN      =  2
ISPIND     =  2
LCHARG     =  True
LHFCALC    =  True
LMAXMIX    =  4
LORBIT     =  11
LPLANE     =  True
LREAL      =  Auto
LSCALU     =  False
LWAVE      =  True
MAGMOM     =  1*6.0 2*-6.0 1*6.0 20*0.6
NKRED      =  2
NPAR       =  8
NSIM       =  1
NSW        =  99
NUPDOWN    =  0
PREC       =  Accurate
SIGMA      =  0.05
SYSTEM     =  Id=[0] dblock_code=[97763-icsd] formula=[li mn (p o4)] sg_name=[p n m a]
TIME       =  0.4"""
        assert s == expected

    def test_lsorbit_magmom(self):
        magmom1 = [[0.0, 0.0, 3.0], [0, 1, 0], [2, 1, 2]]
        magmom2 = [-1, -1, -1, 0, 0, 0, 0, 0]
        magmom4 = [Magmom([1.0, 2.0, 2.0])]

        ans_string1 = "LANGEVIN_GAMMA = 10 10 10\nLSORBIT = True\nMAGMOM = 0.0 0.0 3.0 0 1 0 2 1 2\n"
        ans_string2 = "LANGEVIN_GAMMA = 10\nLSORBIT = True\nMAGMOM = 3*3*-1 3*5*0\n"
        ans_string3 = "LSORBIT = False\nMAGMOM = 2*-1 2*9\n"
        ans_string4_nolsorbit = "LANGEVIN_GAMMA = 10\nLSORBIT = False\nMAGMOM = 1*3.0\n"
        ans_string4_lsorbit = "LANGEVIN_GAMMA = 10\nLSORBIT = True\nMAGMOM = 1.0 2.0 2.0\n"

        incar = Incar({})
        incar["MAGMOM"] = magmom1
        incar["LSORBIT"] = "T"
        incar["LANGEVIN_GAMMA"] = [10, 10, 10]
        assert ans_string1 == str(incar)

        incar["MAGMOM"] = magmom2
        incar["LSORBIT"] = "T"
        incar["LANGEVIN_GAMMA"] = 10
        assert ans_string2 == str(incar)

        incar["MAGMOM"] = magmom4
        incar["LSORBIT"] = "F"
        assert ans_string4_nolsorbit == str(incar)
        incar["LSORBIT"] = "T"
        assert ans_string4_lsorbit == str(incar)

        incar = Incar.from_str(ans_string1)
        assert incar["MAGMOM"] == [[0.0, 0.0, 3.0], [0, 1, 0], [2, 1, 2]]
        assert incar["LANGEVIN_GAMMA"] == [10, 10, 10]

        incar = Incar.from_str(ans_string2)
        assert incar["MAGMOM"] == [
            [-1, -1, -1],
            [-1, -1, -1],
            [-1, -1, -1],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
        ]
        assert incar["LANGEVIN_GAMMA"] == [10]

        incar = Incar.from_str(ans_string3)
        assert not incar["LSORBIT"]
        assert incar["MAGMOM"] == [-1, -1, 9, 9]

    def test_quad_efg(self):
        incar1 = Incar({})
        incar1["LEFG"] = True
        incar1["QUAD_EFG"] = [0.0, 146.6, -25.58]
        ans_string1 = "LEFG = True\nQUAD_EFG = 0.0 146.6 -25.58\n"
        assert ans_string1 == str(incar1)
        incar2 = Incar.from_str(ans_string1)
        assert ans_string1 == str(incar2)

    def test_types(self):
        incar_str = """ALGO = Fast
ECUT = 510
EDIFF = 1e-07
EINT = -0.85 0.85
IBRION = -1
ICHARG = 11
ISIF = 3
ISMEAR = 1
ISPIN = 1
LPARD = True
NBMOD = -3
PREC = Accurate
SIGMA = 0.1"""
        i = Incar.from_str(incar_str)
        assert isinstance(i["EINT"], list)
        assert i["EINT"][0] == -0.85

        incar_str += "\nLHFCALC = .TRUE. ; HFSCREEN = 0.2"
        incar_str += "\nALGO = All;"
        i = Incar.from_str(incar_str)
        assert i["LHFCALC"]
        assert i["HFSCREEN"] == 0.2
        assert i["ALGO"] == "All"

    def test_proc_types(self):
        assert Incar.proc_val("HELLO", "-0.85 0.85") == "-0.85 0.85"

    def test_check_params(self):
        # Triggers warnings when running into nonsensical parameters
        with pytest.warns(BadIncarWarning):
            incar = Incar(
                {
                    "ADDGRID": True,
                    "ALGO": "Normal",
                    "AMIN": 0.01,
                    "AMIX": 0.2,
                    "BMIX": 0.001,
                    "EDIFF": 5 + 1j,  # EDIFF needs to be real
                    "EDIFFG": -0.01,
                    "ENCUT": 520,
                    "IBRION": 2,
                    "ICHARG": 1,
                    "ISIF": 3,
                    "ISMEAR": 1,
                    "ISPIN": 2,
                    "LASPH": 5,  # Should be a bool
                    "LORBIT": 11,
                    "LREAL": "Auto",
                    "LWAVE": False,
                    "MAGMOM": [1, 2, 4, 5],
                    "METAGGA": "SCAM",  # spelling mistake
                    "NELM": 200,
                    "NPAR": 4,
                    "NSW": 99,
                    "PREC": "Accurate",
                    "SIGMA": 0.2,
                    "NBAND": 250,  # spelling mistake
                    "PHON_TLIST": "is_a_str",  # this parameter should be a list
                    "LATTICE_CONSTRAINTS": [
                        True,
                        False,
                        "f",
                    ],  # Should be a list of bools
                    "M_CONSTR": [True, 1, "string"],  # Should be a list of real numbers
                }
            )
            incar.check_params()


class TestKpoints:
    def test_init(self):
        filepath = f"{TEST_FILES_DIR}/KPOINTS.auto"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.kpts == [[10]], "Wrong kpoint lattice read"
        filepath = f"{TEST_FILES_DIR}/KPOINTS.cartesian"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.kpts == [[0.25, 0, 0], [0, 0.25, 0], [0, 0, 0.25]], "Wrong kpoint lattice read"
        assert kpoints.kpts_shift == [0.5, 0.5, 0.5], "Wrong kpoint shift read"

        filepath = f"{TEST_FILES_DIR}/KPOINTS"
        kpoints = Kpoints.from_file(filepath)
        self.kpoints = kpoints
        assert kpoints.kpts == [[2, 4, 6]]

        filepath = f"{TEST_FILES_DIR}/KPOINTS.band"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.labels is not None
        assert kpoints.style == Kpoints.supported_modes.Line_mode
        kpoints_str = str(kpoints)
        assert kpoints_str.split("\n")[3] == "Reciprocal"

        filepath = f"{TEST_FILES_DIR}/KPOINTS.explicit"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.kpts_weights is not None
        assert (
            str(kpoints).strip()
            == """Example file
4
Cartesian
0.0 0.0 0.0 1 None
0.0 0.0 0.5 1 None
0.0 0.5 0.5 2 None
0.5 0.5 0.5 4 None"""
        )

        filepath = f"{TEST_FILES_DIR}/KPOINTS.explicit_tet"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.tet_connections == [(6, [1, 2, 3, 4])]

    def test_style_setter(self):
        filepath = f"{TEST_FILES_DIR}/KPOINTS"
        kpoints = Kpoints.from_file(filepath)
        assert kpoints.style == Kpoints.supported_modes.Monkhorst
        kpoints.style = "G"
        assert kpoints.style == Kpoints.supported_modes.Gamma

    def test_static_constructors(self):
        kpoints = Kpoints.gamma_automatic([3, 3, 3], [0, 0, 0])
        assert kpoints.style == Kpoints.supported_modes.Gamma
        assert kpoints.kpts == [[3, 3, 3]]
        kpoints = Kpoints.monkhorst_automatic([2, 2, 2], [0, 0, 0])
        assert kpoints.style == Kpoints.supported_modes.Monkhorst
        assert kpoints.kpts == [[2, 2, 2]]
        kpoints = Kpoints.automatic(100)
        assert kpoints.style == Kpoints.supported_modes.Automatic
        assert kpoints.kpts == [[100]]
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        poscar = Poscar.from_file(filepath)
        kpoints = Kpoints.automatic_density(poscar.structure, 500)
        assert kpoints.kpts == [[1, 3, 3]]
        assert kpoints.style == Kpoints.supported_modes.Gamma
        kpoints = Kpoints.automatic_density(poscar.structure, 500, True)
        assert kpoints.style == Kpoints.supported_modes.Gamma
        kpoints = Kpoints.automatic_density_by_vol(poscar.structure, 1000)
        assert kpoints.kpts == [[6, 10, 13]]
        assert kpoints.style == Kpoints.supported_modes.Gamma
        kpoints = Kpoints.automatic_density_by_lengths(poscar.structure, [50, 50, 1], True)
        assert kpoints.kpts == [[5, 9, 1]]
        assert kpoints.style == Kpoints.supported_modes.Gamma

        struct = poscar.structure
        struct.make_supercell(3)
        kpoints = Kpoints.automatic_density(struct, 500)
        assert kpoints.kpts == [[1, 1, 1]]
        assert kpoints.style == Kpoints.supported_modes.Gamma
        kpoints = Kpoints.from_str(
            """k-point mesh
0
G
10 10 10
0.5 0.5 0.5
"""
        )
        assert np.allclose(kpoints.kpts_shift, [0.5, 0.5, 0.5])

    def test_as_dict_from_dict(self):
        k = Kpoints.monkhorst_automatic([2, 2, 2], [0, 0, 0])
        d = k.as_dict()
        k2 = Kpoints.from_dict(d)
        assert k.kpts == k2.kpts
        assert k.style == k2.style
        assert k.kpts_shift == k2.kpts_shift

    def test_kpt_bands_as_dict_from_dict(self):
        file_name = f"{TEST_FILES_DIR}/KPOINTS.band"
        k = Kpoints.from_file(file_name)
        d = k.as_dict()
        import json

        json.dumps(d)
        # This doesn't work
        k2 = Kpoints.from_dict(d)
        assert k.kpts == k2.kpts
        assert k.style == k2.style
        assert k.kpts_shift == k2.kpts_shift
        assert k.num_kpts == k2.num_kpts

    def test_pickle(self):
        k = Kpoints.gamma_automatic()
        pickle.dumps(k)

    def test_automatic_kpoint(self):
        # struct = PymatgenTest.get_structure("Li2O")
        p = Poscar.from_str(
            """Al1
1.0
2.473329 0.000000 1.427977
0.824443 2.331877 1.427977
0.000000 0.000000 2.855955
Al
1
direct
0.000000 0.000000 0.000000 Al"""
        )
        kpoints = Kpoints.automatic_density(p.structure, 1000)
        assert np.allclose(kpoints.kpts[0], [10, 10, 10])

    def test_automatic_density_by_lengths(self):
        # Load a structure from a POSCAR file
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        structure = Poscar.from_file(filepath).structure

        # test different combos of length densities and expected kpoints
        # TODO should test Monkhorst style case and force_gamma=True case
        for length_densities, expected_kpts, expected_style in [
            ([50, 50, 1], [[5, 9, 1]], Kpoints.supported_modes.Gamma),
            ([25, 50, 3], [[3, 9, 1]], Kpoints.supported_modes.Gamma),
            ([24, 48, 2], [[3, 8, 1]], Kpoints.supported_modes.Gamma),
        ]:
            kpoints = Kpoints.automatic_density_by_lengths(structure, length_densities)

            assert kpoints.kpts == expected_kpts

            assert kpoints.style == expected_style


class TestPotcarSingle:
    _multiprocess_shared_ = True

    def setup(self):
        self.psingle = PotcarSingle.from_file(f"{TEST_FILES_DIR}/POT_GGA_PAW_PBE/POTCAR.Mn_pv.gz")

    def test_keywords(self):
        data = {
            "VRHFIN": "Mn: 3p4s3d",
            "LPAW": True,
            "DEXC": -0.003,
            "STEP": [20.000, 1.050],
            "RPACOR": 2.080,
            "LEXCH": "PE",
            "ENMAX": 269.865,
            "QCUT": -4.454,
            "TITEL": "PAW_PBE Mn_pv 07Sep2000",
            "LCOR": True,
            "EAUG": 569.085,
            "RMAX": 2.807,
            "ZVAL": 13.000,
            "EATOM": 2024.8347,
            "NDATA": 100,
            "LULTRA": False,
            "QGAM": 8.907,
            "ENMIN": 202.399,
            "RCLOC": 1.725,
            "RCORE": 2.300,
            "RDEP": 2.338,
            "IUNSCR": 1,
            "RAUG": 1.300,
            "POMASS": 54.938,
            "RWIGS": 1.323,
        }
        assert self.psingle.keywords == data

    def test_psctr(self):
        filename = f"{TEST_FILES_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe.gz"

        psingle = PotcarSingle.from_file(filename)

        data = {
            "nentries": 9,
            "Orbitals": (
                (1, 0, 0.50, -6993.8440, 2.0000),
                (2, 0, 0.50, -0814.6047, 2.0000),
                (2, 1, 1.50, -0693.3689, 6.0000),
                (3, 0, 0.50, -0089.4732, 2.0000),
                (3, 1, 1.50, -0055.6373, 6.0000),
                (3, 2, 2.50, -0003.8151, 7.0000),
                (4, 0, 0.50, -0004.2551, 1.0000),
                (4, 1, 1.50, -0003.4015, 0.0000),
                (4, 3, 2.50, -0001.3606, 0.0000),
            ),
            "OrbitalDescriptions": (
                (2, -3.8151135, 23, 2.300, None, None),
                (2, -5.1756961, 23, 2.300, None, None),
                (0, -4.2550963, 23, 2.300, None, None),
                (0, 07.2035603, 23, 2.300, None, None),
                (1, -2.7211652, 23, 2.300, None, None),
                (1, 18.4316424, 23, 2.300, None, None),
            ),
        }
        for k, v in data.items():
            assert psingle.PSCTR[k] == v

    def test_nelectrons(self):
        assert self.psingle.nelectrons == 13

    def test_electron_config(self):
        config = self.psingle.electron_configuration
        assert config[-1] == (3, "p", 6)

    def test_attributes(self):
        for k in [
            "DEXC",
            "RPACOR",
            "ENMAX",
            "QCUT",
            "EAUG",
            "RMAX",
            "ZVAL",
            "EATOM",
            "NDATA",
            "QGAM",
            "ENMIN",
            "RCLOC",
            "RCORE",
            "RDEP",
            "RAUG",
            "POMASS",
            "RWIGS",
        ]:
            assert getattr(self.psingle, k) is not None

    def test_found_unknown_key(self):
        with pytest.raises(KeyError, match="BAD_KEY"):
            PotcarSingle.parse_functions["BAD_KEY"]

    def test_bad_value(self):
        with pytest.raises(ValueError, match="could not convert string to float"):
            PotcarSingle.parse_functions["ENMAX"]("ThisShouldBeAFloat")

    def test_hash(self):
        assert self.psingle.get_potcar_hash() == "fa52f891f234d49bb4cb5ea96aae8f98"

    def test_functional_types(self):
        assert self.psingle.functional == "PBE"

        assert self.psingle.functional_class == "GGA"

        assert self.psingle.potential_type == "PAW"

        psingle = PotcarSingle.from_file(f"{TEST_FILES_DIR}/POT_LDA_PAW/POTCAR.Fe.gz")

        assert psingle.functional == "Perdew-Zunger81"

        assert psingle.functional_class == "LDA"

        assert psingle.potential_type == "PAW"

        assert self.psingle.symbol == "Mn_pv"

    def test_identify_potcar(self):
        filename = f"{TEST_FILES_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe.gz"

        psingle = PotcarSingle.from_file(filename)
        assert "PBE_54" in psingle.identify_potcar()[0]
        assert "Fe" in psingle.identify_potcar()[1]

    def test_potcar_hash_warning(self):
        filename = f"{TEST_FILES_DIR}/modified_potcars_data/POT_GGA_PAW_PBE/POTCAR.Fe_pv"
        with pytest.warns(UnknownPotcarWarning, match="POTCAR is known to match the following functionals:"):
            PotcarSingle.from_file(filename)

    def test_potcar_file_hash_warning(self):
        filename = f"{TEST_FILES_DIR}/modified_potcars_header/POT_GGA_PAW_PBE/POTCAR.Fe_pv"
        with pytest.warns(UnknownPotcarWarning, match="POTCAR is corrupted"):
            PotcarSingle.from_file(filename)

    def test_verify_faulty_potcar_with_hash(self):
        filename = f"{TEST_FILES_DIR}/modified_potcars_data/POT_GGA_PAW_PBE_54/POTCAR.Fe_pv_with_hash"
        with pytest.warns(UnknownPotcarWarning, match="POTCAR with symbol Fe_pv has metadata that "):
            PotcarSingle.from_file(filename)

    def test_verify_correct_potcar_with_hash(self):
        filename = f"{TEST_FILES_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe_pv_with_hash.gz"
        vaspdir = os.path.abspath(os.path.dirname(pymatgen.io.vasp.__file__))
        file_hash_db = loadfn(f"{vaspdir}/vasp_potcar_file_hashes.json")
        metadata_hash_db = loadfn(f"{vaspdir}/vasp_potcar_pymatgen_hashes.json")

        psingle = PotcarSingle.from_file(filename)
        assert psingle.hash in metadata_hash_db
        assert psingle.file_hash in file_hash_db
        assert psingle.hash_sha256_computed == psingle.hash_sha256_from_file

    def test_multi_potcar_with_and_without_hash(self):
        filename = f"{TEST_FILES_DIR}/POT_GGA_PAW_PBE_54/POTCAR.Fe_O.gz"
        vaspdir = os.path.abspath(os.path.dirname(pymatgen.io.vasp.__file__))
        loadfn(f"{vaspdir}/vasp_potcar_file_hashes.json")
        Potcar.from_file(filename)
        # Still need to test the if POTCAR can be read.
        # No longer testing for hashes
        # for psingle in potcars:
        #     if hasattr(psingle, "hash_sha256_from_file"):
        #         assert psingle.hash_sha256_computed == psingle.hash_sha256_from_file
        # else:
        #     assert psingle.file_hash in file_hash_db

    # def test_default_functional(self):
    #     p = PotcarSingle.from_symbol_and_functional("Fe")
    #     assert p.functional_class == "GGA"
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "LDA"
    #     p = PotcarSingle.from_symbol_and_functional("Fe")
    #     assert p.functional_class == "LDA"
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "PBE"


class TestPotcar:
    def setup(self):
        if "PMG_VASP_PSP_DIR" not in SETTINGS:
            SETTINGS["PMG_VASP_PSP_DIR"] = str(TEST_FILES_DIR)
        self.filepath = f"{TEST_FILES_DIR}/POTCAR"
        self.potcar = Potcar.from_file(self.filepath)

    def test_init(self):
        assert self.potcar.symbols == ["Fe", "P", "O"], "Wrong symbols read in for POTCAR"
        potcar = Potcar(["Fe_pv", "O"])
        assert potcar[0].enmax == 293.238

    def test_from_file(self):
        assert {d.header for d in self.potcar} == {"PAW_PBE O 08Apr2002", "PAW_PBE P 17Jan2003", "PAW_PBE Fe 06Sep2000"}

    def test_potcar_map(self):
        fe_potcar = zopen(f"{TEST_FILES_DIR}/POT_GGA_PAW_PBE/POTCAR.Fe_pv.gz").read().decode("utf-8")
        # specify V instead of Fe - this makes sure the test won't pass if the
        # code just grabs the POTCAR from the config file (the config file would
        # grab the V POTCAR)
        potcar = Potcar(["V"], sym_potcar_map={"V": fe_potcar})
        assert potcar.symbols == ["Fe_pv"], "Wrong symbols read in for POTCAR"

    def test_to_from_dict(self):
        d = self.potcar.as_dict()
        potcar = Potcar.from_dict(d)
        assert potcar.symbols == ["Fe", "P", "O"]

    def test_write(self):
        tempfname = Path("POTCAR.testing")
        self.potcar.write_file(tempfname)
        p = Potcar.from_file(tempfname)
        assert p.symbols == self.potcar.symbols

        # check line by line
        with open(self.filepath) as f_ref, open(tempfname) as f_new:
            ref_potcar = f_ref.readlines()
            new_potcar = f_new.readlines()

        if len(ref_potcar) != len(new_potcar):
            raise AssertionError("POTCAR file lengths are not equal")

        for line1, line2 in zip(ref_potcar, new_potcar):
            if line1.strip() != line2.strip():
                raise AssertionError("POTCAR contents are not")

        tempfname.unlink()

    def test_set_symbol(self):
        assert self.potcar.symbols == ["Fe", "P", "O"]
        assert self.potcar[0].nelectrons == 8
        self.potcar.symbols = ["Fe_pv", "O"]
        assert self.potcar.symbols == ["Fe_pv", "O"]
        assert self.potcar[0].nelectrons == 14

    # def test_default_functional(self):
    #     p = Potcar(["Fe", "P"])
    #     assert p[0].functional_class == "GGA"
    #     assert p[1].functional_class == "GGA"
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "LDA"
    #     p = Potcar(["Fe", "P"])
    #     assert p[0].functional_class == "LDA"
    #     assert p[1].functional_class == "LDA"

    def test_pickle(self):
        pickle.dumps(self.potcar)

    # def tearDown(self):
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "PBE"


class TestVaspInput:
    def setup(self):
        filepath = f"{TEST_FILES_DIR}/INCAR"
        incar = Incar.from_file(filepath)
        filepath = f"{TEST_FILES_DIR}/POSCAR"
        poscar = Poscar.from_file(filepath, check_for_POTCAR=False)
        if "PMG_VASP_PSP_DIR" not in os.environ:
            os.environ["PMG_VASP_PSP_DIR"] = str(TEST_FILES_DIR)
        filepath = f"{TEST_FILES_DIR}/POTCAR"
        potcar = Potcar.from_file(filepath)
        filepath = f"{TEST_FILES_DIR}/KPOINTS.auto"
        kpoints = Kpoints.from_file(filepath)
        self.vasp_input = VaspInput(incar, kpoints, poscar, potcar)

    def test_to_from_dict(self):
        d = self.vasp_input.as_dict()
        vasp_input = VaspInput.from_dict(d)
        comp = vasp_input["POSCAR"].structure.composition
        assert comp == Composition("Fe4P4O16")

    def test_write(self):
        tmp_dir = Path("VaspInput.testing")
        self.vasp_input.write_input(tmp_dir)

        filepath = tmp_dir / "INCAR"
        incar = Incar.from_file(filepath)
        assert incar["NSW"] == 99

        for name in ("INCAR", "POSCAR", "POTCAR", "KPOINTS"):
            (tmp_dir / name).unlink()

        tmp_dir.rmdir()

    def test_run_vasp(self):
        self.vasp_input.run_vasp(".", vasp_cmd=["cat", "INCAR"])
        with open("vasp.out") as f:
            output = f.read()
            assert output.split("\n")[0] == "ALGO = Damped"

    def test_from_directory(self):
        vi = VaspInput.from_directory(TEST_FILES_DIR, optional_files={"CONTCAR.Li2O": Poscar})
        assert vi["INCAR"]["ALGO"] == "Damped"
        assert "CONTCAR.Li2O" in vi
        d = vi.as_dict()
        vasp_input = VaspInput.from_dict(d)
        assert "CONTCAR.Li2O" in vasp_input
