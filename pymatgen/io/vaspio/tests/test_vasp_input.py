#!/usr/bin/env python

'''
Created on Jul 16, 2012
'''

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Jul 16, 2012"

import unittest
import os
import numpy as np

from pymatgen.core.physical_constants import AMU_TO_KG, BOLTZMANN_CONST
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Kpoints, Potcar, PotcarSingle
from pymatgen import Composition, Structure, __file__

test_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'test_files')

class PoscarTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        comp = poscar.structure.composition
        self.assertEqual(comp, Composition.from_formula("Fe4P4O16"))

        #Vasp 4 type with symbols at the end.
        poscar_string = """Test1
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
direct
0.000000 0.000000 0.000000 Si
0.750000 0.500000 0.750000 F"""
        poscar = Poscar.from_string(poscar_string)
        self.assertEqual(poscar.structure.composition, Composition.from_formula("SiF"))

        #Vasp 4 tyle file with default names, i.e. no element symbol found.
        poscar_string = """Test2
1.0
3.840198 0.000000 0.000000
1.920099 3.325710 0.000000
0.000000 -2.217138 3.135509
1 1
direct
0.000000 0.000000 0.000000
0.750000 0.500000 0.750000"""
        poscar = Poscar.from_string(poscar_string)
        self.assertEqual(poscar.structure.composition, Composition.from_formula("HHe"))

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
0.750000 0.500000 0.750000 F F F O"""
        poscar = Poscar.from_string(poscar_string)
        self.assertEqual(poscar.selective_dynamics, [[True, True, True], [False, False, False]])
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
0.750000 0.500000 0.750000 F F F O"""
        poscar = Poscar.from_string(poscar_string)
        d = poscar.to_dict
        poscar2 = Poscar.from_dict(d)
        self.assertEqual(poscar2.comment, "Test3")
        self.assertTrue(all(poscar2.selective_dynamics[0]))
        self.assertFalse(all(poscar2.selective_dynamics[1]))

    def test_str(self):
        si = 14
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])

        #Silicon structure for testing.
        latt = [[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]]
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
0.750000 0.500000 0.750000 Si'''

        self.assertEquals(str(poscar), expected_str, "Wrong POSCAR output!")

    def test_from_md_run(self):
        #Parsing from an MD type run with velocities
        p = Poscar.from_file(os.path.join(test_dir, "CONTCAR.MD"), check_for_POTCAR=False)
        self.assertAlmostEqual(np.sum(np.array(p.velocities)), 0.0065417961324)
        self.assertEqual(p.predictor_corrector[0][0], 1)
        self.assertEqual(p.predictor_corrector[1][0], 2)

    def test_setattr(self):
        filepath = os.path.join(test_dir, 'POSCAR')
        poscar = Poscar.from_file(filepath)
        self.assertRaises(ValueError, setattr, poscar, 'velocities', [[0, 0, 0]])

    def test_velocities(self):
        si = 14
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.75, 0.5, 0.75])

        #Silicon structure for testing.
        latt = [[ 3.8401979337, 0.00, 0.00], [1.9200989668, 3.3257101909, 0.00], [0.00, -2.2171384943, 3.1355090603]]
        struct = Structure(latt, [si, si], coords)
        poscar = Poscar(struct)
        poscar.set_temperature(900)

        v = np.array(poscar.velocities)

        for x in np.sum(v, axis=0):
            self.assertAlmostEqual(x, 0, 7, 'Velocities initialized with a net momentum')

        temperature = struct[0].specie.atomic_mass * AMU_TO_KG * np.sum(v ** 2) / (3 * BOLTZMANN_CONST) * 1e10
        self.assertAlmostEqual(temperature, 900, 4, 'Temperature instantiated incorrectly')


class IncarTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'INCAR')
        incar = Incar.from_file(filepath)
        incar["LDAU"] = "T"
        self.assertEqual(incar["ALGO"], "Damped", "Wrong Algo")
        self.assertEqual(float(incar["EDIFF"]), 1e-4, "Wrong EDIFF")

    def test_diff(self):
        filepath1 = os.path.join(test_dir, 'INCAR')
        incar1 = Incar.from_file(filepath1)
        filepath2 = os.path.join(test_dir, 'INCAR.2')
        incar2 = Incar.from_file(filepath2)
        self.assertEqual(incar1.diff(incar2), {'Different': {'NELM': {'INCAR1': 'Default', 'INCAR2': 100}, 'ISPIND': {'INCAR1': 2, 'INCAR2': 'Default'}, 'LWAVE': {'INCAR1': True, 'INCAR2': False}, 'LDAUPRINT': {'INCAR1': 'Default', 'INCAR2': 1}, 'MAGMOM': {'INCAR1': [6, -6, -6, 6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6], 'INCAR2': 'Default'}, 'NELMIN': {'INCAR1': 'Default', 'INCAR2': 3}, 'ENCUTFOCK': {'INCAR1': 0.0, 'INCAR2': 'Default'}, 'HFSCREEN': {'INCAR1': 0.207, 'INCAR2': 'Default'}, 'LSCALU': {'INCAR1': False, 'INCAR2': 'Default'}, 'ENCUT': {'INCAR1': 500, 'INCAR2': 'Default'}, 'NSIM': {'INCAR1': 1, 'INCAR2': 'Default'}, 'ICHARG': {'INCAR1': 'Default', 'INCAR2': 1}, 'NSW': {'INCAR1': 99, 'INCAR2': 51}, 'NKRED': {'INCAR1': 2, 'INCAR2': 'Default'}, 'NUPDOWN': {'INCAR1': 0, 'INCAR2': 'Default'}, 'LCHARG': {'INCAR1': True, 'INCAR2': 'Default'}, 'LPLANE': {'INCAR1': True, 'INCAR2': 'Default'}, 'ISMEAR': {'INCAR1': 0, 'INCAR2':-5}, 'NPAR': {'INCAR1': 8, 'INCAR2': 1}, 'SYSTEM': {'INCAR1': 'Id=[0] dblock_code=[97763-icsd] formula=[li mn (p o4)] sg_name=[p n m a]', 'INCAR2': 'Id=[91090] dblock_code=[20070929235612linio-59.53134651-vasp] formula=[li3 ni3 o6] sg_name=[r-3m]'}, 'ALGO': {'INCAR1': 'Damped', 'INCAR2': 'Fast'}, 'LHFCALC': {'INCAR1': True, 'INCAR2': 'Default'}, 'TIME': {'INCAR1': 0.4, 'INCAR2': 'Default'}}, 'Same': {'IBRION': 2, 'PREC': 'Accurate', 'ISIF': 3, 'LMAXMIX': 4, 'LREAL': 'Auto', 'ISPIN': 2, 'EDIFF': 0.0001, 'LORBIT': '11', 'SIGMA': 0.05}})

    def test_to_dict_and_from_dict(self):
        file_name = os.path.join(test_dir, 'INCAR')
        incar = Incar.from_file(file_name)
        d = incar.to_dict
        incar2 = Incar.from_dict(d)
        self.assertEqual(incar, incar2)


class KpointsTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'KPOINTS.auto')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts, [[10]], "Wrong kpoint lattice read")
        filepath = os.path.join(test_dir, 'KPOINTS.cartesian')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts, [[0.25, 0, 0], [0, 0.25, 0], [0, 0, 0.25]], "Wrong kpoint lattice read")
        self.assertEqual(kpoints.kpts_shift, [0.5, 0.5, 0.5], "Wrong kpoint shift read")

        filepath = os.path.join(test_dir, 'KPOINTS')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts, [[2, 4, 6]], "Wrong kpoint lattice read")

        filepath = os.path.join(test_dir, 'KPOINTS.band')
        kpoints = Kpoints.from_file(filepath)
        self.assertIsNotNone(kpoints.labels)
        self.assertEqual(kpoints.style, "Line_mode")

        filepath = os.path.join(test_dir, 'KPOINTS.explicit')
        kpoints = Kpoints.from_file(filepath)
        self.assertIsNotNone(kpoints.kpts_weights)

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

    def test_to_dict_from_dict(self):
        k = Kpoints.monkhorst_automatic([2, 2, 2], [0, 0, 0])
        d = k.to_dict
        k2 = Kpoints.from_dict(d)
        self.assertEqual(k.kpts, k2.kpts)
        self.assertEqual(k.style, k2.style)
        self.assertEqual(k.kpts_shift, k2.kpts_shift)

    def test_kpt_bands_to_dict_from_dict(self):
        file_name = os.path.join(test_dir, 'KPOINTS.band')
        k = Kpoints.from_file(file_name)
        d = k.to_dict
        #This doesn't work
        #k2 = Kpoints.from_dict(d)
        #self.assertEqual(k.kpts, k2.kpts)
        #self.assertEqual(k.style, k2.style)
        #self.assertEqual(k.kpts_shift, k2.kpts_shift)
        #self.assertEqual(k.num_kpts, k2.num_kpts)


class PotcarSingleTest(unittest.TestCase):

    def setUp(self):
        with open(os.path.join(test_dir, 'POTCAR.Mn_pv'), 'r') as f:
            self.psingle = PotcarSingle(f.read())

    def test_keywords(self):
        data = {'VRHFIN': 'Mn: 3p4s3d', 'LPAW': 'T    paw PP', 'DEXC': '-.003',
                'STEP': '20.000   1.050',
                'RPACOR': '2.080    partial core radius', 'LEXCH': 'PE',
                'ENMAX': '269.865', 'QCUT': '-4.454',
                'TITEL': 'PAW_PBE Mn_pv 07Sep2000',
                'LCOR': 'T    correct aug charges', 'EAUG': '569.085',
                'RMAX': '2.807    core radius for proj-oper',
                'ZVAL': '13.000    mass and valenz',
                'EATOM': '2024.8347 eV,  148.8212 Ry', 'NDATA': '100',
                'LULTRA': 'F    use ultrasoft PP ?',
                'QGAM': '8.907    optimization parameters',
                'ENMIN': '202.399 eV', 'RCLOC': '1.725    cutoff for local pot',
                'RCORE': '2.300    outmost cutoff radius',
                'RDEP': '2.338    radius for radial grids',
                'IUNSCR': '1    unscreen: 0-lin 1-nonlin 2-no',
                'RAUG': '1.300    factor for augmentation sphere',
                'POMASS': '54.938',
                'RWIGS': '1.323    wigner-seitz radius (au A)'}
        self.assertEqual(self.psingle.keywords, data)

    def test_nelectrons(self):
        self.assertEqual(self.psingle.nelectrons, 13)

    def test_attributes(self):
        for k in ['DEXC', 'RPACOR', 'ENMAX', 'QCUT', 'EAUG', 'RMAX',
                         'ZVAL', 'EATOM', 'NDATA', 'QGAM', 'ENMIN', 'RCLOC',
                         'RCORE', 'RDEP', 'RAUG', 'POMASS', 'RWIGS']:
            self.assertIsNotNone(getattr(self.psingle, k))

    def test_from_functional_and_symbols(self):
        p = PotcarSingle.from_symbol_and_functional("Li_sv", "PBE")
        self.assertEqual(p.enmax, 271.649)

class PotcarTest(unittest.TestCase):

    def test_init(self):
        filepath = os.path.join(test_dir, 'POTCAR')
        potcar = Potcar.from_file(filepath)
        self.assertEqual(potcar.symbols, ["Fe", "P", "O"], "Wrong symbols read in for POTCAR")
        potcar = Potcar(["Fe_pv", "O"])
        self.assertEqual(potcar[0].enmax, 293.238)

    def test_potcar_map(self):
        fe_potcar = open(os.path.join(test_dir, 'Fe_POTCAR')).read()
        #specify V instead of Fe - this makes sure the test won't pass if the
        #code just grabs the POTCAR from the config file (the config file would
        #grab the V POTCAR)
        potcar = Potcar(["V"], sym_potcar_map={"V": fe_potcar})
        self.assertEqual(potcar.symbols, ["Fe"], "Wrong symbols read in for POTCAR")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
