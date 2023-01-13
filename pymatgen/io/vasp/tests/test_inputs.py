# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import annotations

import os
import pickle
import unittest
import warnings
from pathlib import Path

import numpy as np
import pytest
import scipy.constants as const
from monty.io import zopen
from monty.serialization import loadfn
from monty.tempfile import ScratchDir

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
from pymatgen.util.testing import PymatgenTest


class PoscarTest(PymatgenTest):
    def test_init(self):
        filepath = PymatgenTest.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath, check_for_POTCAR=False)
        comp = poscar.structure.composition
        self.assertEqual(comp, Composition("Fe4P4O16"))

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
        poscar = Poscar.from_string(poscar_string)
        self.assertEqual(poscar.structure.composition, Composition("SiF"))

        poscar_string = ""
        self.assertRaises(ValueError, Poscar.from_string, poscar_string)

        # VASP 4 tyle file with default names, i.e. no element symbol found.
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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            poscar = Poscar.from_string(poscar_string)
        self.assertEqual(poscar.structure.composition, Composition("HHe"))
        # VASP 4 tyle file with default names, i.e. no element symbol found.
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
        poscar = Poscar.from_string(poscar_string)
        self.assertEqual(poscar.selective_dynamics, [[True, True, True], [False, False, False]])
        self.selective_poscar = poscar

    def test_from_file(self):
        filepath = PymatgenTest.TEST_FILES_DIR / "POSCAR.symbols_natoms_multilines"
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
        self.assertEqual([site.specie.symbol for site in poscar.structure], ordered_expected_elements)

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
        poscar = Poscar.from_string(poscar_string)
        d = poscar.as_dict()
        poscar2 = Poscar.from_dict(d)
        self.assertEqual(poscar2.comment, "Test3")
        self.assertTrue(all(poscar2.selective_dynamics[0]))
        self.assertFalse(all(poscar2.selective_dynamics[1]))

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
        p = Poscar.from_string(poscar_string)
        site = p.structure[1]
        self.assertArrayAlmostEqual(site.coords, np.array([3.840198, 1.5, 2.35163175]) * 1.1)

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
        self.assertEqual(actual_str, expected_str, "Wrong POSCAR output!")

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

        self.assertEqual(str(poscar), expected_str, "Wrong POSCAR output!")

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
        poscar = Poscar.from_string(poscar_string)
        self.assertEqual(str(poscar), expected)

    def test_from_md_run(self):
        # Parsing from an MD type run with velocities and predictor corrector data
        p = Poscar.from_file(PymatgenTest.TEST_FILES_DIR / "CONTCAR.MD", check_for_POTCAR=False)
        self.assertAlmostEqual(np.sum(np.array(p.velocities)), 0.0065417961324)
        self.assertEqual(p.predictor_corrector[0][0][0], 0.33387820e00)
        self.assertEqual(p.predictor_corrector[0][1][1], -0.10583589e-02)

    def test_write_MD_poscar(self):
        # Parsing from an MD type run with velocities and predictor corrector data
        # And writing a new POSCAR from the new structure
        p = Poscar.from_file(PymatgenTest.TEST_FILES_DIR / "CONTCAR.MD", check_for_POTCAR=False)

        tempfname = Path("POSCAR.testing.md")
        p.write_file(tempfname)
        p3 = Poscar.from_file(tempfname)

        self.assertArrayAlmostEqual(p.structure.lattice.abc, p3.structure.lattice.abc, 5)
        self.assertArrayAlmostEqual(p.velocities, p3.velocities, 5)
        self.assertArrayAlmostEqual(p.predictor_corrector, p3.predictor_corrector, 5)
        self.assertEqual(p.predictor_corrector_preamble, p3.predictor_corrector_preamble)
        tempfname.unlink()

    def test_setattr(self):
        filepath = PymatgenTest.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath, check_for_POTCAR=False)
        self.assertRaises(ValueError, setattr, poscar, "velocities", [[0, 0, 0]])
        poscar.selective_dynamics = np.array([[True, False, False]] * 24)
        ans = """
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
        self.assertEqual(str(poscar).strip(), ans.strip())

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
            self.assertAlmostEqual(x, 0, 7)

        temperature = struct[0].specie.atomic_mass.to("kg") * np.sum(v**2) / (3 * const.k) * 1e10
        self.assertAlmostEqual(temperature, 900, 4, "Temperature instantiated incorrectly")

        poscar.set_temperature(700)
        v = np.array(poscar.velocities)
        for x in np.sum(v, axis=0):
            self.assertAlmostEqual(x, 0, 7, "Velocities initialized with a net momentum")

        temperature = struct[0].specie.atomic_mass.to("kg") * np.sum(v**2) / (3 * const.k) * 1e10
        self.assertAlmostEqual(temperature, 700, 4, "Temperature instantiated incorrectly")

    def test_write(self):
        filepath = PymatgenTest.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath)
        tempfname = Path("POSCAR.testing")
        poscar.write_file(tempfname)
        p = Poscar.from_file(tempfname)
        self.assertArrayAlmostEqual(poscar.structure.lattice.abc, p.structure.lattice.abc, 5)
        tempfname.unlink()


class IncarTest(PymatgenTest):
    def setUp(self):
        file_name = PymatgenTest.TEST_FILES_DIR / "INCAR"
        self.incar = Incar.from_file(file_name)

    def test_init(self):
        incar = self.incar
        incar["LDAU"] = "T"
        self.assertEqual(incar["ALGO"], "Damped", "Wrong Algo")
        self.assertEqual(float(incar["EDIFF"]), 1e-4, "Wrong EDIFF")
        self.assertEqual(type(incar["LORBIT"]), int)

    def test_diff(self):
        filepath1 = PymatgenTest.TEST_FILES_DIR / "INCAR"
        incar1 = Incar.from_file(filepath1)
        filepath2 = PymatgenTest.TEST_FILES_DIR / "INCAR.2"
        incar2 = Incar.from_file(filepath2)
        incar3 = Incar.from_file(filepath2)
        self.assertEqual(
            incar1.diff(incar2),
            {
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
            },
        )

        self.assertEqual(
            incar1.diff(incar3),
            {
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
            },
        )

    def test_as_dict_and_from_dict(self):
        d = self.incar.as_dict()
        incar2 = Incar.from_dict(d)
        self.assertEqual(self.incar, incar2)
        d["MAGMOM"] = [Magmom([1, 2, 3]).as_dict()]
        incar3 = Incar.from_dict(d)
        self.assertEqual(incar3["MAGMOM"], [Magmom([1, 2, 3])])

    def test_write(self):
        tempfname = Path("INCAR.testing")
        self.incar.write_file(tempfname)
        i = Incar.from_file(tempfname)
        self.assertEqual(i, self.incar)
        tempfname.unlink()

    def test_get_string(self):
        s = self.incar.get_string(pretty=True, sort_keys=True)
        ans = """ALGO       =  Damped
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
        self.assertEqual(s, ans)

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
        self.assertEqual(ans_string1, str(incar))

        incar["MAGMOM"] = magmom2
        incar["LSORBIT"] = "T"
        incar["LANGEVIN_GAMMA"] = 10
        self.assertEqual(ans_string2, str(incar))

        incar["MAGMOM"] = magmom4
        incar["LSORBIT"] = "F"
        self.assertEqual(ans_string4_nolsorbit, str(incar))
        incar["LSORBIT"] = "T"
        self.assertEqual(ans_string4_lsorbit, str(incar))

        incar = Incar.from_string(ans_string1)
        self.assertEqual(incar["MAGMOM"], [[0.0, 0.0, 3.0], [0, 1, 0], [2, 1, 2]])
        self.assertEqual(incar["LANGEVIN_GAMMA"], [10, 10, 10])

        incar = Incar.from_string(ans_string2)
        self.assertEqual(
            incar["MAGMOM"],
            [
                [-1, -1, -1],
                [-1, -1, -1],
                [-1, -1, -1],
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
            ],
        )
        self.assertEqual(incar["LANGEVIN_GAMMA"], [10])

        incar = Incar.from_string(ans_string3)
        self.assertFalse(incar["LSORBIT"])
        self.assertEqual(incar["MAGMOM"], [-1, -1, 9, 9])

    def test_quad_efg(self):
        incar1 = Incar({})
        incar1["LEFG"] = True
        incar1["QUAD_EFG"] = [0.0, 146.6, -25.58]
        ans_string1 = "LEFG = True\nQUAD_EFG = 0.0 146.6 -25.58\n"
        self.assertEqual(ans_string1, str(incar1))
        incar2 = Incar.from_string(ans_string1)
        self.assertEqual(ans_string1, str(incar2))

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
        i = Incar.from_string(incar_str)
        self.assertIsInstance(i["EINT"], list)
        self.assertEqual(i["EINT"][0], -0.85)

        incar_str += "\nLHFCALC = .TRUE. ; HFSCREEN = 0.2"
        incar_str += "\nALGO = All;"
        i = Incar.from_string(incar_str)
        self.assertTrue(i["LHFCALC"])
        self.assertEqual(i["HFSCREEN"], 0.2)
        self.assertEqual(i["ALGO"], "All")

    def test_proc_types(self):
        self.assertEqual(Incar.proc_val("HELLO", "-0.85 0.85"), "-0.85 0.85")

    def test_check_params(self):
        # Triggers warnings when running into nonsensical parameters
        with self.assertWarns(BadIncarWarning):
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


class KpointsTest(PymatgenTest):
    def test_init(self):
        filepath = PymatgenTest.TEST_FILES_DIR / "KPOINTS.auto"
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts, [[10]], "Wrong kpoint lattice read")
        filepath = PymatgenTest.TEST_FILES_DIR / "KPOINTS.cartesian"
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(
            kpoints.kpts,
            [[0.25, 0, 0], [0, 0.25, 0], [0, 0, 0.25]],
            "Wrong kpoint lattice read",
        )
        self.assertEqual(kpoints.kpts_shift, [0.5, 0.5, 0.5], "Wrong kpoint shift read")

        filepath = PymatgenTest.TEST_FILES_DIR / "KPOINTS"
        kpoints = Kpoints.from_file(filepath)
        self.kpoints = kpoints
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])

        filepath = PymatgenTest.TEST_FILES_DIR / "KPOINTS.band"
        kpoints = Kpoints.from_file(filepath)
        self.assertIsNotNone(kpoints.labels)
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Line_mode)
        kpoints_str = str(kpoints)
        self.assertEqual(kpoints_str.split("\n")[3], "Reciprocal")

        filepath = PymatgenTest.TEST_FILES_DIR / "KPOINTS.explicit"
        kpoints = Kpoints.from_file(filepath)
        self.assertIsNotNone(kpoints.kpts_weights)
        self.assertEqual(
            str(kpoints).strip(),
            """Example file
4
Cartesian
0.0 0.0 0.0 1 None
0.0 0.0 0.5 1 None
0.0 0.5 0.5 2 None
0.5 0.5 0.5 4 None""",
        )

        filepath = PymatgenTest.TEST_FILES_DIR / "KPOINTS.explicit_tet"
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.tet_connections, [(6, [1, 2, 3, 4])])

    def test_style_setter(self):
        filepath = PymatgenTest.TEST_FILES_DIR / "KPOINTS"
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)
        kpoints.style = "G"
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

    def test_static_constructors(self):
        kpoints = Kpoints.gamma_automatic([3, 3, 3], [0, 0, 0])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)
        self.assertEqual(kpoints.kpts, [[3, 3, 3]])
        kpoints = Kpoints.monkhorst_automatic([2, 2, 2], [0, 0, 0])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Monkhorst)
        self.assertEqual(kpoints.kpts, [[2, 2, 2]])
        kpoints = Kpoints.automatic(100)
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Automatic)
        self.assertEqual(kpoints.kpts, [[100]])
        filepath = PymatgenTest.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath)
        kpoints = Kpoints.automatic_density(poscar.structure, 500)
        self.assertEqual(kpoints.kpts, [[1, 3, 3]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)
        kpoints = Kpoints.automatic_density(poscar.structure, 500, True)
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)
        kpoints = Kpoints.automatic_density_by_vol(poscar.structure, 1000)
        self.assertEqual(kpoints.kpts, [[6, 10, 13]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)
        kpoints = Kpoints.automatic_density_by_lengths(poscar.structure, [50, 50, 1], True)
        self.assertEqual(kpoints.kpts, [[5, 9, 1]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)

        s = poscar.structure
        s.make_supercell(3)
        kpoints = Kpoints.automatic_density(s, 500)
        self.assertEqual(kpoints.kpts, [[1, 1, 1]])
        self.assertEqual(kpoints.style, Kpoints.supported_modes.Gamma)
        kpoints = Kpoints.from_string(
            """k-point mesh
0
G
10 10 10
0.5 0.5 0.5
"""
        )
        self.assertArrayAlmostEqual(kpoints.kpts_shift, [0.5, 0.5, 0.5])

    def test_as_dict_from_dict(self):
        k = Kpoints.monkhorst_automatic([2, 2, 2], [0, 0, 0])
        d = k.as_dict()
        k2 = Kpoints.from_dict(d)
        self.assertEqual(k.kpts, k2.kpts)
        self.assertEqual(k.style, k2.style)
        self.assertEqual(k.kpts_shift, k2.kpts_shift)

    def test_kpt_bands_as_dict_from_dict(self):
        file_name = PymatgenTest.TEST_FILES_DIR / "KPOINTS.band"
        k = Kpoints.from_file(file_name)
        d = k.as_dict()
        import json

        json.dumps(d)
        # This doesn't work
        k2 = Kpoints.from_dict(d)
        self.assertEqual(k.kpts, k2.kpts)
        self.assertEqual(k.style, k2.style)
        self.assertEqual(k.kpts_shift, k2.kpts_shift)
        self.assertEqual(k.num_kpts, k2.num_kpts)

    def test_pickle(self):
        k = Kpoints.gamma_automatic()
        pickle.dumps(k)

    def test_automatic_kpoint(self):
        # s = PymatgenTest.get_structure("Li2O")
        p = Poscar.from_string(
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
        self.assertArrayAlmostEqual(kpoints.kpts[0], [10, 10, 10])


class PotcarSingleTest(PymatgenTest):
    _multiprocess_shared_ = True

    def setUp(self):
        self.psingle = PotcarSingle.from_file(PymatgenTest.TEST_FILES_DIR / "POT_GGA_PAW_PBE" / "POTCAR.Mn_pv.gz")

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
        self.assertEqual(self.psingle.keywords, data)

    def test_psctr(self):
        filename = PymatgenTest.TEST_FILES_DIR / "POT_GGA_PAW_PBE_54" / "POTCAR.Fe.gz"

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
            self.assertEqual(psingle.PSCTR[k], v)

    def test_nelectrons(self):
        self.assertEqual(self.psingle.nelectrons, 13)

    def test_electron_config(self):
        config = self.psingle.electron_configuration
        self.assertEqual(config[-1], (3, "p", 6))

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
            self.assertIsNotNone(getattr(self.psingle, k))

    def test_found_unknown_key(self):
        with self.assertRaises(KeyError):
            PotcarSingle.parse_functions["BAD_KEY"]

    def test_bad_value(self):
        self.assertRaises(ValueError, PotcarSingle.parse_functions["ENMAX"], "ThisShouldBeAFloat")

    def test_hash(self):
        self.assertEqual(self.psingle.get_potcar_hash(), "fa52f891f234d49bb4cb5ea96aae8f98")

    def test_functional_types(self):
        self.assertEqual(self.psingle.functional, "PBE")

        self.assertEqual(self.psingle.functional_class, "GGA")

        self.assertEqual(self.psingle.potential_type, "PAW")

        psingle = PotcarSingle.from_file(PymatgenTest.TEST_FILES_DIR / "POT_LDA_PAW" / "POTCAR.Fe.gz")

        self.assertEqual(psingle.functional, "Perdew-Zunger81")

        self.assertEqual(psingle.functional_class, "LDA")

        self.assertEqual(psingle.potential_type, "PAW")

    def test_identify_potcar(self):
        filename = PymatgenTest.TEST_FILES_DIR / "POT_GGA_PAW_PBE_54" / "POTCAR.Fe.gz"

        psingle = PotcarSingle.from_file(filename)
        assert "PBE_54" in psingle.identify_potcar()[0]
        assert "Fe" in psingle.identify_potcar()[1]

    def test_potcar_hash_warning(self):
        filename = PymatgenTest.TEST_FILES_DIR / "modified_potcars_data" / "POT_GGA_PAW_PBE" / "POTCAR.Fe_pv"
        with pytest.warns(UnknownPotcarWarning, match="POTCAR is known to match the following functionals:"):
            PotcarSingle.from_file(filename)

    def test_potcar_file_hash_warning(self):
        filename = PymatgenTest.TEST_FILES_DIR / "modified_potcars_header" / "POT_GGA_PAW_PBE" / "POTCAR.Fe_pv"
        with pytest.warns(UnknownPotcarWarning, match="POTCAR is corrupted"):
            PotcarSingle.from_file(filename)

    def test_verify_faulty_potcar_with_hash(self):
        filename = (
            PymatgenTest.TEST_FILES_DIR / "modified_potcars_data" / "POT_GGA_PAW_PBE_54" / "POTCAR.Fe_pv_with_hash"
        )
        with pytest.raises(ValueError):
            PotcarSingle.from_file(filename)

    def test_verify_correct_potcar_with_hash(self):
        filename = PymatgenTest.TEST_FILES_DIR / "POT_GGA_PAW_PBE_54" / "POTCAR.Fe_pv_with_hash.gz"
        cwd = os.path.abspath(os.path.dirname(__file__))
        file_hash_db = loadfn(os.path.join(cwd, "../vasp_potcar_file_hashes.json"))
        metadata_hash_db = loadfn(os.path.join(cwd, "../vasp_potcar_pymatgen_hashes.json"))

        psingle = PotcarSingle.from_file(filename)
        assert psingle.hash in metadata_hash_db
        assert psingle.file_hash in file_hash_db
        assert psingle.hash_sha256_computed == psingle.hash_sha256_from_file

    def test_multi_potcar_with_and_without_hash(self):
        filename = PymatgenTest.TEST_FILES_DIR / "POT_GGA_PAW_PBE_54" / "POTCAR.Fe_O.gz"
        cwd = os.path.abspath(os.path.dirname(__file__))
        file_hash_db = loadfn(os.path.join(cwd, "../vasp_potcar_file_hashes.json"))
        potcars = Potcar.from_file(filename)
        for psingle in potcars:
            if hasattr(psingle, "hash_sha256_from_file"):
                assert psingle.hash_sha256_computed == psingle.hash_sha256_from_file
            else:
                assert psingle.file_hash in file_hash_db

    # def test_default_functional(self):
    #     p = PotcarSingle.from_symbol_and_functional("Fe")
    #     self.assertEqual(p.functional_class, 'GGA')
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "LDA"
    #     p = PotcarSingle.from_symbol_and_functional("Fe")
    #     self.assertEqual(p.functional_class, 'LDA')
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "PBE"


class PotcarTest(PymatgenTest):
    def setUp(self):
        if "PMG_VASP_PSP_DIR" not in SETTINGS:
            SETTINGS["PMG_VASP_PSP_DIR"] = str(PymatgenTest.TEST_FILES_DIR)
        filepath = PymatgenTest.TEST_FILES_DIR / "POTCAR"
        self.potcar = Potcar.from_file(filepath)

    def test_init(self):
        self.assertEqual(self.potcar.symbols, ["Fe", "P", "O"], "Wrong symbols read in for POTCAR")
        potcar = Potcar(["Fe_pv", "O"])
        self.assertEqual(potcar[0].enmax, 293.238)

    def test_potcar_map(self):
        fe_potcar = zopen(PymatgenTest.TEST_FILES_DIR / "POT_GGA_PAW_PBE" / "POTCAR.Fe_pv.gz").read().decode("utf-8")
        # specify V instead of Fe - this makes sure the test won't pass if the
        # code just grabs the POTCAR from the config file (the config file would
        # grab the V POTCAR)
        potcar = Potcar(["V"], sym_potcar_map={"V": fe_potcar})
        self.assertEqual(potcar.symbols, ["Fe_pv"], "Wrong symbols read in for POTCAR")

    def test_to_from_dict(self):
        d = self.potcar.as_dict()
        potcar = Potcar.from_dict(d)
        self.assertEqual(potcar.symbols, ["Fe", "P", "O"])

    def test_write(self):
        tempfname = Path("POTCAR.testing")
        self.potcar.write_file(tempfname)
        p = Potcar.from_file(tempfname)
        self.assertEqual(p.symbols, self.potcar.symbols)
        tempfname.unlink()

    def test_set_symbol(self):
        self.assertEqual(self.potcar.symbols, ["Fe", "P", "O"])
        self.assertEqual(self.potcar[0].nelectrons, 8)
        self.potcar.symbols = ["Fe_pv", "O"]
        self.assertEqual(self.potcar.symbols, ["Fe_pv", "O"])
        self.assertEqual(self.potcar[0].nelectrons, 14)

    # def test_default_functional(self):
    #     p = Potcar(["Fe", "P"])
    #     self.assertEqual(p[0].functional_class, 'GGA')
    #     self.assertEqual(p[1].functional_class, 'GGA')
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "LDA"
    #     p = Potcar(["Fe", "P"])
    #     self.assertEqual(p[0].functional_class, 'LDA')
    #     self.assertEqual(p[1].functional_class, 'LDA')

    def test_pickle(self):
        pickle.dumps(self.potcar)

    # def tearDown(self):
    #     SETTINGS["PMG_DEFAULT_FUNCTIONAL"] = "PBE"


class VaspInputTest(PymatgenTest):
    def setUp(self):
        filepath = PymatgenTest.TEST_FILES_DIR / "INCAR"
        incar = Incar.from_file(filepath)
        filepath = PymatgenTest.TEST_FILES_DIR / "POSCAR"
        poscar = Poscar.from_file(filepath, check_for_POTCAR=False)
        if "PMG_VASP_PSP_DIR" not in os.environ:
            os.environ["PMG_VASP_PSP_DIR"] = str(PymatgenTest.TEST_FILES_DIR)
        filepath = PymatgenTest.TEST_FILES_DIR / "POTCAR"
        potcar = Potcar.from_file(filepath)
        filepath = PymatgenTest.TEST_FILES_DIR / "KPOINTS.auto"
        kpoints = Kpoints.from_file(filepath)
        self.vinput = VaspInput(incar, kpoints, poscar, potcar)

    def test_to_from_dict(self):
        d = self.vinput.as_dict()
        vinput = VaspInput.from_dict(d)
        comp = vinput["POSCAR"].structure.composition
        self.assertEqual(comp, Composition("Fe4P4O16"))

    def test_write(self):
        tmp_dir = Path("VaspInput.testing")
        self.vinput.write_input(tmp_dir)

        filepath = tmp_dir / "INCAR"
        incar = Incar.from_file(filepath)
        self.assertEqual(incar["NSW"], 99)

        for name in ("INCAR", "POSCAR", "POTCAR", "KPOINTS"):
            (tmp_dir / name).unlink()

        tmp_dir.rmdir()

    def test_run_vasp(self):
        # To add some test.
        with ScratchDir(".") as d:
            self.vinput.run_vasp(d, vasp_cmd=["cat", "INCAR"])
            with open(os.path.join(d, "vasp.out")) as f:
                output = f.read()
                self.assertEqual(output.split("\n")[0], "ALGO = Damped")

    def test_from_directory(self):
        vi = VaspInput.from_directory(PymatgenTest.TEST_FILES_DIR, optional_files={"CONTCAR.Li2O": Poscar})
        self.assertEqual(vi["INCAR"]["ALGO"], "Damped")
        self.assertIn("CONTCAR.Li2O", vi)
        d = vi.as_dict()
        vinput = VaspInput.from_dict(d)
        self.assertIn("CONTCAR.Li2O", vinput)


if __name__ == "__main__":
    unittest.main()
