# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Jul 16, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 16, 2012"

import unittest
import os
import numpy as np
import warnings

from pymatgen.util.testing import PymatgenTest
from pymatgen.core.physical_constants import BOLTZMANN_CONST
from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints, Potcar, \
    PotcarSingle, VaspInput
from pymatgen import Composition, Structure
from monty.io import zopen


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                        'test_files')


class PoscarTest(PymatgenTest):

    def test_init(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        comp = poscar.structure.composition
        self.assertEqual(comp, Composition("Fe4P4O16"))

        #Vasp 4 type with symbols at the end.
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

        #Vasp 4 tyle file with default names, i.e. no element symbol found.
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
        #Vasp 4 tyle file with default names, i.e. no element symbol found.
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
        self.assertEqual(poscar.selective_dynamics, [[True, True, True],
                                                     [False, False, False]])
        self.selective_poscar = poscar

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
        self.assertArrayAlmostEqual(site.coords,
                                    np.array([3.840198, 1.5, 2.35163175]) * 1.1)

    def test_str(self):
        si = 14
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])

        #Silicon structure for testing.
        latt = [[3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(latt, [si, si], coords)
        poscar = Poscar(struct)
        expected_str = '''Si2
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
Si
2
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 Si
'''

        self.assertEqual(str(poscar), expected_str, "Wrong POSCAR output!")

        #Vasp 4 type with symbols at the end.
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
3.840198 -0.000000 -0.000000
-1.920099 -3.325710 -0.000000
-0.000000 2.217138 -3.135509
Si F
1 1
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 F
"""
        poscar = Poscar.from_string(poscar_string)
        self.assertEqual(str(poscar), expected)

    def test_from_md_run(self):
        #Parsing from an MD type run with velocities
        p = Poscar.from_file(os.path.join(test_dir, "CONTCAR.MD"),
                             check_for_POTCAR=False)
        self.assertAlmostEqual(np.sum(np.array(p.velocities)), 0.0065417961324)
        self.assertEqual(p.predictor_corrector[0][0], 1)
        self.assertEqual(p.predictor_corrector[1][0], 2)

    def test_setattr(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        self.assertRaises(ValueError, setattr, poscar, 'velocities',
                          [[0, 0, 0]])
        poscar.selective_dynamics = np.array([[True, False, False]] * 24)

    def test_velocities(self):
        si = 14
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])

        #Silicon structure for testing.
        latt = [[3.8401979337, 0.00, 0.00],
                [1.9200989668, 3.3257101909, 0.00],
                [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(latt, [si, si], coords)
        poscar = Poscar(struct)
        poscar.set_temperature(900)

        v = np.array(poscar.velocities)

        for x in np.sum(v, axis=0):
            self.assertAlmostEqual(
                x, 0, 7, 'Velocities initialized with a net momentum')

        temperature = struct[0].specie.atomic_mass.to("kg") * \
            np.sum(v ** 2) / (3 * BOLTZMANN_CONST) * 1e10
        self.assertAlmostEqual(temperature, 900, 4,
                               'Temperature instantiated incorrectly')


    def test_write(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        tempfname = "POSCAR.testing"
        poscar.write_file(tempfname)
        p = Poscar.from_file(tempfname)
        self.assertArrayAlmostEqual(poscar.structure.lattice.abc,
                                    p.structure.lattice.abc, 5)
        os.remove(tempfname)

class IncarTest(unittest.TestCase):

    def setUp(self):
        file_name = os.path.join(test_dir, 'INCAR')
        self.incar = Incar.from_file(file_name)


    def test_init(self):
        incar = self.incar
        incar["LDAU"] = "T"
        self.assertEqual(incar["ALGO"], "Damped", "Wrong Algo")
        self.assertEqual(float(incar["EDIFF"]), 1e-4, "Wrong EDIFF")
        self.assertEqual(type(incar["LORBIT"]), int)

    def test_diff(self):
        incar = self.incar
        filepath1 = os.path.join(test_dir, 'INCAR')
        incar1 = Incar.from_file(filepath1)
        filepath2 = os.path.join(test_dir, 'INCAR.2')
        incar2 = Incar.from_file(filepath2)
        filepath3 = os.path.join(test_dir, 'INCAR.3')
        incar3 = Incar.from_file(filepath2)
        self.assertEqual(
            incar1.diff(incar2),
            {'Different': {
                'NELM': {'INCAR1': None, 'INCAR2': 100},
                'ISPIND': {'INCAR1': 2, 'INCAR2': None},
                'LWAVE': {'INCAR1': True, 'INCAR2': False},
                'LDAUPRINT': {'INCAR1': None, 'INCAR2': 1},
                'MAGMOM': {'INCAR1': [6, -6, -6, 6, 0.6, 0.6, 0.6,
                                      0.6, 0.6, 0.6, 0.6, 0.6,
                                      0.6, 0.6, 0.6, 0.6, 0.6,
                                      0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6],
                           'INCAR2': None},
                'NELMIN': {'INCAR1': None, 'INCAR2': 3},
                'ENCUTFOCK': {'INCAR1': 0.0, 'INCAR2': None},
                'HFSCREEN': {'INCAR1': 0.207, 'INCAR2': None},
                'LSCALU': {'INCAR1': False, 'INCAR2': None},
                'ENCUT': {'INCAR1': 500, 'INCAR2': None},
                'NSIM': {'INCAR1': 1, 'INCAR2': None},
                'ICHARG': {'INCAR1': None, 'INCAR2': 1},
                'NSW': {'INCAR1': 99, 'INCAR2': 51},
                'NKRED': {'INCAR1': 2, 'INCAR2': None},
                'NUPDOWN': {'INCAR1': 0, 'INCAR2': None},
                'LCHARG': {'INCAR1': True, 'INCAR2': None},
                'LPLANE': {'INCAR1': True, 'INCAR2': None},
                'ISMEAR': {'INCAR1': 0, 'INCAR2': -5},
                'NPAR': {'INCAR1': 8, 'INCAR2': 1},
                'SYSTEM': {'INCAR1': 'Id=[0] dblock_code=[97763-icsd] formula=[li mn (p o4)] sg_name=[p n m a]',
                           'INCAR2': 'Id=[91090] dblock_code=[20070929235612linio-59.53134651-vasp] formula=[li3 ni3 o6] sg_name=[r-3m]'},
                'ALGO': {'INCAR1': 'Damped', 'INCAR2': 'Fast'},
                'LHFCALC': {'INCAR1': True, 'INCAR2': None},
                'TIME': {'INCAR1': 0.4, 'INCAR2': None}},
             'Same': {'IBRION': 2, 'PREC': 'Accurate', 'ISIF': 3, 'LMAXMIX': 4,
                      'LREAL': 'Auto', 'ISPIN': 2, 'EDIFF': 0.0001,
                      'LORBIT': 11, 'SIGMA': 0.05}})

        self.assertEqual(
            incar1.diff(incar3),
            {'Different': {
                'NELM': {'INCAR1': None, 'INCAR2': 100},
                'ISPIND': {'INCAR1': 2, 'INCAR2': None},
                'LWAVE': {'INCAR1': True, 'INCAR2': False},
                'LDAUPRINT': {'INCAR1': None, 'INCAR2': 1},
                'MAGMOM': {'INCAR1': [6, -6, -6, 6, 0.6, 0.6, 0.6,
                                      0.6, 0.6, 0.6, 0.6, 0.6,
                                      0.6, 0.6, 0.6, 0.6, 0.6,
                                      0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6],
                           'INCAR2': None},
                'NELMIN': {'INCAR1': None, 'INCAR2': 3},
                'ENCUTFOCK': {'INCAR1': 0.0, 'INCAR2': None},
                'HFSCREEN': {'INCAR1': 0.207, 'INCAR2': None},
                'LSCALU': {'INCAR1': False, 'INCAR2': None},
                'ENCUT': {'INCAR1': 500, 'INCAR2': None},
                'NSIM': {'INCAR1': 1, 'INCAR2': None},
                'ICHARG': {'INCAR1': None, 'INCAR2': 1},
                'NSW': {'INCAR1': 99, 'INCAR2': 51},
                'NKRED': {'INCAR1': 2, 'INCAR2': None},
                'NUPDOWN': {'INCAR1': 0, 'INCAR2': None},
                'LCHARG': {'INCAR1': True, 'INCAR2': None},
                'LPLANE': {'INCAR1': True, 'INCAR2': None},
                'ISMEAR': {'INCAR1': 0, 'INCAR2': -5},
                'NPAR': {'INCAR1': 8, 'INCAR2': 1},
                'SYSTEM': {'INCAR1': 'Id=[0] dblock_code=[97763-icsd] formula=[li mn (p o4)] sg_name=[p n m a]',
                           'INCAR2': 'Id=[91090] dblock_code=[20070929235612linio-59.53134651-vasp] formula=[li3 ni3 o6] sg_name=[r-3m]'},
                'ALGO': {'INCAR1': 'Damped', 'INCAR2': 'Fast'},
                'LHFCALC': {'INCAR1': True, 'INCAR2': None},
                'TIME': {'INCAR1': 0.4, 'INCAR2': None}},
             'Same': {'IBRION': 2, 'PREC': 'Accurate', 'ISIF': 3, 'LMAXMIX': 4,
                      'LREAL': 'Auto', 'ISPIN': 2, 'EDIFF': 0.0001,
                      'LORBIT': 11, 'SIGMA': 0.05}})

    def test_as_dict_and_from_dict(self):
        d = self.incar.as_dict()
        incar2 = Incar.from_dict(d)
        self.assertEqual(self.incar, incar2)

    def test_write(self):
        tempfname = "INCAR.testing"
        self.incar.write_file(tempfname)
        i = Incar.from_file(tempfname)
        self.assertEqual(i, self.incar)
        os.remove(tempfname)


class KpointsTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'KPOINTS.auto')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts, [[10]], "Wrong kpoint lattice read")
        filepath = os.path.join(test_dir, 'KPOINTS.cartesian')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts,
                         [[0.25, 0, 0], [0, 0.25, 0], [0, 0, 0.25]],
                         "Wrong kpoint lattice read")
        self.assertEqual(kpoints.kpts_shift, [0.5, 0.5, 0.5],
                         "Wrong kpoint shift read")

        filepath = os.path.join(test_dir, 'KPOINTS')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts, [[2, 4, 6]])

        filepath = os.path.join(test_dir, 'KPOINTS.band')
        kpoints = Kpoints.from_file(filepath)
        self.assertIsNotNone(kpoints.labels)
        self.assertEqual(kpoints.style, "Line_mode")
        kpoints_str = str(kpoints)
        self.assertEqual(kpoints_str.split("\n")[3], "Reciprocal")

        filepath = os.path.join(test_dir, 'KPOINTS.explicit')
        kpoints = Kpoints.from_file(filepath)
        self.assertIsNotNone(kpoints.kpts_weights)
        self.assertEqual(str(kpoints).strip(), """Example file
4
Cartesian
0.0 0.0 0.0 1 None
0.0 0.0 0.5 1 None
0.0 0.5 0.5 2 None
0.5 0.5 0.5 4 None""")

        filepath = os.path.join(test_dir, 'KPOINTS.explicit_tet')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.tet_connections, [(6, [1, 2, 3, 4])])

    def test_static_constructors(self):
        kpoints = Kpoints.gamma_automatic([3, 3, 3], [0, 0, 0])
        self.assertEqual(kpoints.style, "Gamma")
        self.assertEqual(kpoints.kpts, [[3, 3, 3]])
        kpoints = Kpoints.monkhorst_automatic([2, 2, 2], [0, 0, 0])
        self.assertEqual(kpoints.style, "Monkhorst")
        self.assertEqual(kpoints.kpts, [[2, 2, 2]])
        kpoints = Kpoints.automatic(100)
        self.assertEqual(kpoints.style, "Automatic")
        self.assertEqual(kpoints.kpts, [[100]])
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        kpoints = Kpoints.automatic_density(poscar.structure, 500)
        self.assertEqual(kpoints.kpts, [[2, 4, 4]])
        self.assertEqual(kpoints.style, "Monkhorst")
        kpoints = Kpoints.automatic_density(poscar.structure, 500, True)
        self.assertEqual(kpoints.style, "Gamma")
        kpoints = Kpoints.automatic_density_by_vol(poscar.structure, 1000)
        self.assertEqual(kpoints.kpts, [[6, 11, 13]])
        self.assertEqual(kpoints.style, "Gamma")

        s = poscar.structure
        s.make_supercell(3)
        kpoints = Kpoints.automatic_density(s, 500)
        self.assertEqual(kpoints.kpts, [[1, 1, 1]])
        self.assertEqual(kpoints.style, "Gamma")

    def test_as_dict_from_dict(self):
        k = Kpoints.monkhorst_automatic([2, 2, 2], [0, 0, 0])
        d = k.as_dict()
        k2 = Kpoints.from_dict(d)
        self.assertEqual(k.kpts, k2.kpts)
        self.assertEqual(k.style, k2.style)
        self.assertEqual(k.kpts_shift, k2.kpts_shift)

    def test_kpt_bands_as_dict_from_dict(self):
        file_name = os.path.join(test_dir, 'KPOINTS.band')
        k = Kpoints.from_file(file_name)
        d = k.as_dict()
        import json
        json.dumps(d)
        #This doesn't work
        k2 = Kpoints.from_dict(d)
        self.assertEqual(k.kpts, k2.kpts)
        self.assertEqual(k.style, k2.style)
        self.assertEqual(k.kpts_shift, k2.kpts_shift)
        self.assertEqual(k.num_kpts, k2.num_kpts)


class PotcarSingleTest(unittest.TestCase):

    def setUp(self):
        #with zopen(os.path.join(test_dir, "POT_GGA_PAW_PBE",
        #                        "POTCAR.Mn_pv.gz"), 'rb') as f:
        self.psingle = PotcarSingle.from_file(os.path.join(test_dir, "POT_GGA_PAW_PBE",
                                "POTCAR.Mn_pv.gz"))

    def test_keywords(self):
        data = {'VRHFIN': 'Mn: 3p4s3d', 'LPAW': True, 'DEXC': -.003,
                'STEP': [20.000,   1.050],
                'RPACOR': 2.080, 'LEXCH': 'PE',
                'ENMAX': 269.865, 'QCUT': -4.454,
                'TITEL': 'PAW_PBE Mn_pv 07Sep2000',
                'LCOR': True, 'EAUG': 569.085,
                'RMAX': 2.807,
                'ZVAL': 13.000,
                'EATOM': 2024.8347, 'NDATA': 100,
                'LULTRA': False,
                'QGAM': 8.907,
                'ENMIN': 202.399,
                'RCLOC': 1.725,
                'RCORE': 2.300,
                'RDEP': 2.338,
                'IUNSCR': 1,
                'RAUG': 1.300,
                'POMASS': 54.938,
                'RWIGS': 1.323}
        self.assertEqual(self.psingle.keywords, data)

    def test_nelectrons(self):
        self.assertEqual(self.psingle.nelectrons, 13)

    def test_electron_config(self):
        config = self.psingle.electron_configuration
        self.assertEqual(config[-1], (3, "p", 6))

    def test_attributes(self):
        for k in ['DEXC', 'RPACOR', 'ENMAX', 'QCUT', 'EAUG', 'RMAX',
                  'ZVAL', 'EATOM', 'NDATA', 'QGAM', 'ENMIN', 'RCLOC',
                  'RCORE', 'RDEP', 'RAUG', 'POMASS', 'RWIGS']:
            self.assertIsNotNone(getattr(self.psingle, k))

    def test_found_unknown_key(self):
        with self.assertRaises(KeyError):
            PotcarSingle.parse_functions['BAD_KEY']

    def test_bad_value(self):
        self.assertRaises(ValueError, PotcarSingle.parse_functions['ENMAX'], "ThisShouldBeAFloat")

    def test_hash(self):
        self.assertEqual(self.psingle.get_potcar_hash(), "fa52f891f234d49bb4cb5ea96aae8f98")

    def test_from_functional_and_symbols(self):
        if "VASP_PSP_DIR" not in os.environ:
            test_potcar_dir = os.path.abspath(
                os.path.join(os.path.dirname(__file__),
                             "..", "..", "..", "..", "test_files"))
            os.environ["VASP_PSP_DIR"] = test_potcar_dir
        p = PotcarSingle.from_symbol_and_functional("Li_sv", "PBE")
        self.assertEqual(p.enmax, 271.649)

    def test_functional_types(self):
        self.assertEqual(self.psingle.functional, 'PBE')

        self.assertEqual(self.psingle.functional_class, 'GGA')

        self.assertEqual(self.psingle.potential_type, 'PAW')

        psingle = PotcarSingle.from_file(os.path.join(test_dir, "POT_LDA_PAW",
                                "POTCAR.Fe.gz"))

        self.assertEqual(psingle.functional, 'Perdew-Zunger81')

        self.assertEqual(psingle.functional_class, 'LDA')

        self.assertEqual(psingle.potential_type, 'PAW')


class PotcarTest(unittest.TestCase):

    def setUp(self):
        if "VASP_PSP_DIR" not in os.environ:
            test_potcar_dir = os.path.abspath(
                os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                             "test_files"))
            os.environ["VASP_PSP_DIR"] = test_potcar_dir
        filepath = os.path.join(test_dir, 'POTCAR')
        self.potcar = Potcar.from_file(filepath)

    def test_init(self):
        self.assertEqual(self.potcar.symbols, ["Fe", "P", "O"],
                         "Wrong symbols read in for POTCAR")
        potcar = Potcar(["Fe_pv", "O"])
        self.assertEqual(potcar[0].enmax, 293.238)

    def test_potcar_map(self):
        fe_potcar = zopen(os.path.join(test_dir, "POT_GGA_PAW_PBE",
                                       "POTCAR.Fe_pv.gz")).read().decode(
            "utf-8")
        #specify V instead of Fe - this makes sure the test won't pass if the
        #code just grabs the POTCAR from the config file (the config file would
        #grab the V POTCAR)
        potcar = Potcar(["V"], sym_potcar_map={"V": fe_potcar})
        self.assertEqual(potcar.symbols, ["Fe_pv"], "Wrong symbols read in "
                                                    "for POTCAR")

    def test_to_from_dict(self):
        d = self.potcar.as_dict()
        potcar = Potcar.from_dict(d)
        self.assertEqual(potcar.symbols, ["Fe", "P", "O"])

    def test_write(self):
        tempfname = "POTCAR.testing"
        self.potcar.write_file(tempfname)
        p = Potcar.from_file(tempfname)
        self.assertEqual(p.symbols, self.potcar.symbols)
        os.remove(tempfname)

    def test_set_symbol(self):
        self.assertEqual(self.potcar.symbols, ["Fe", "P", "O"])
        self.assertEqual(self.potcar[0].nelectrons, 8)
        self.potcar.symbols = ["Fe_pv", "O"]
        self.assertEqual(self.potcar.symbols, ["Fe_pv", "O"])
        self.assertEqual(self.potcar[0].nelectrons, 14)


class VaspInputTest(unittest.TestCase):

    def setUp(self):
        filepath = os.path.join(test_dir, 'INCAR')
        incar = Incar.from_file(filepath)
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        if "VASP_PSP_DIR" not in os.environ:
            test_potcar_dir = os.path.abspath(
                os.path.join(os.path.dirname(__file__), "..", "..", "..", "..",
                             "test_files"))
            os.environ["VASP_PSP_DIR"] = test_potcar_dir
        filepath = os.path.join(test_dir, 'POTCAR')
        potcar = Potcar.from_file(filepath)
        filepath = os.path.join(test_dir, 'KPOINTS.auto')
        kpoints = Kpoints.from_file(filepath)
        self.vinput = VaspInput(incar, kpoints, poscar, potcar)

    def test_to_from_dict(self):
        d = self.vinput.as_dict()
        vinput = VaspInput.from_dict(d)
        comp = vinput["POSCAR"].structure.composition
        self.assertEqual(comp, Composition("Fe4P4O16"))

    def test_write(self):
        tmp_dir = "VaspInput.testing"
        self.vinput.write_input(tmp_dir)

        filepath = os.path.join(tmp_dir, "INCAR")
        incar = Incar.from_file(filepath)
        self.assertEqual(incar["NSW"], 99)

        for name in ("INCAR", "POSCAR", "POTCAR", "KPOINTS"):
            os.remove(os.path.join(tmp_dir, name))

        os.rmdir(tmp_dir)

    def test_from_directory(self):
        vi = VaspInput.from_directory(test_dir,
                                      optional_files={"CONTCAR.Li2O": Poscar})
        self.assertEqual(vi["INCAR"]["ALGO"], "Damped")
        self.assertIn("CONTCAR.Li2O", vi)
        d = vi.as_dict()
        vinput = VaspInput.from_dict(d)
        self.assertIn("CONTCAR.Li2O", vinput)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
