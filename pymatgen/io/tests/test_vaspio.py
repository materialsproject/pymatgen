#!/usr/bin/python
import unittest
import os

from pymatgen.io.vaspio import Poscar, Potcar, Kpoints, Incar, Vasprun, Outcar, Oszicar
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Composition, Structure
from numpy import array

module_dir = os.path.dirname(os.path.abspath(__file__))

class  PoscarTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','POSCAR')
        poscar = Poscar.from_file(filepath)
        comp = poscar.struct.composition
        self.assertEqual(comp,Composition.from_formula("Fe4P4O16"))
        
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
        self.assertEqual(poscar.struct.composition, Composition.from_formula("SiF"))
        
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
        self.assertEqual(poscar.struct.composition, Composition.from_formula("HHe"))
        
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
        
    def test_str(self):
        si = 14
        coords = list()
        coords.append(array([0,0,0]))
        coords.append(array([0.75,0.5,0.75]))

        #Silicon structure for testing.
        latt = Lattice(array([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]]))
        struct = Structure(latt,[si,si],coords)
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

class  IncarTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','INCAR')
        incar = Incar.from_file(filepath)
        incar["LDAU"] = "T"
        self.assertEqual(incar["ALGO"],"Damped","Wrong Algo")
        self.assertEqual(float(incar["EDIFF"]),1e-4,"Wrong EDIFF")        
        
    def test_diff(self):
        filepath1 = os.path.join(module_dir, 'vasp_testfiles','INCAR')
        incar1 = Incar.from_file(filepath1)
        filepath2 = os.path.join(module_dir, 'vasp_testfiles','INCAR.2')
        incar2 = Incar.from_file(filepath2)
        self.assertEqual(incar1.diff(incar2), {'Different': {'NELM': {'INCAR1': 'Default', 'INCAR2': 100}, 'ISPIND': {'INCAR1': 2, 'INCAR2': 'Default'}, 'LWAVE': {'INCAR1': True, 'INCAR2': False}, 'LDAUPRINT': {'INCAR1': 'Default', 'INCAR2': 1}, 'MAGMOM': {'INCAR1': [6, -6, -6, 6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6], 'INCAR2': 'Default'}, 'NELMIN': {'INCAR1': 'Default', 'INCAR2': 3}, 'ENCUTFOCK': {'INCAR1': 0.0, 'INCAR2': 'Default'}, 'HFSCREEN': {'INCAR1': 0.207, 'INCAR2': 'Default'}, 'LSCALU': {'INCAR1': False, 'INCAR2': 'Default'}, 'ENCUT': {'INCAR1': 500, 'INCAR2': 'Default'}, 'NSIM': {'INCAR1': 1, 'INCAR2': 'Default'}, 'ICHARG': {'INCAR1': 'Default', 'INCAR2': 1}, 'NSW': {'INCAR1': 99, 'INCAR2': 51}, 'NKRED': {'INCAR1': 2, 'INCAR2': 'Default'}, 'NUPDOWN': {'INCAR1': 0, 'INCAR2': 'Default'}, 'LCHARG': {'INCAR1': True, 'INCAR2': 'Default'}, 'LPLANE': {'INCAR1': True, 'INCAR2': 'Default'}, 'ISMEAR': {'INCAR1': 0, 'INCAR2': -5}, 'NPAR': {'INCAR1': 8, 'INCAR2': 1}, 'SYSTEM': {'INCAR1': 'Id=[0] dblock_code=[97763-icsd] formula=[li mn (p o4)] sg_name=[p n m a]', 'INCAR2': 'Id=[91090] dblock_code=[20070929235612linio-59.53134651-vasp] formula=[li3 ni3 o6] sg_name=[r-3m]'}, 'ALGO': {'INCAR1': 'Damped', 'INCAR2': 'Fast'}, 'LHFCALC': {'INCAR1': True, 'INCAR2': 'Default'}, 'TIME': {'INCAR1': 0.4, 'INCAR2': 'Default'}}, 'Same': {'IBRION': 2, 'PREC': 'Accurate', 'ISIF': 3, 'LMAXMIX': 4, 'LREAL': 'Auto', 'ISPIN': 2, 'EDIFF': 0.0001, 'LORBIT': '11', 'SIGMA': 0.05}})
    
class  KpointsTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','KPOINTS.auto')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts,[[10]],"Wrong kpoint lattice read")
        filepath = os.path.join(module_dir, 'vasp_testfiles','KPOINTS.cartesian')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts,[[0.25, 0, 0], [0,0.25,0], [0,0,0.25]],"Wrong kpoint lattice read")
        self.assertEqual(kpoints.kpts_shift,[0.5, 0.5, 0.5],"Wrong kpoint shift read")
        
        filepath = os.path.join(module_dir, 'vasp_testfiles','KPOINTS')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts,[[2,4,6]],"Wrong kpoint lattice read")
        
        filepath = os.path.join(module_dir, 'vasp_testfiles','KPOINTS.band')
        kpoints = Kpoints.from_file(filepath)
        self.assertIsNotNone(kpoints.labels)
        self.assertEqual(kpoints.style, "Line_mode")
        
        filepath = os.path.join(module_dir, 'vasp_testfiles','KPOINTS.explicit')
        kpoints = Kpoints.from_file(filepath)
        self.assertIsNotNone(kpoints.kpts_weights)
        
        filepath = os.path.join(module_dir, 'vasp_testfiles','KPOINTS.explicit_tet')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.tet_connections, [(6, [1,2,3,4])])
        
    def test_static_constructors(self):
        kpoints = Kpoints.gamma_automatic([3,3,3], [0,0,0])
        self.assertEqual(kpoints.style, "Gamma")
        self.assertEqual(kpoints.kpts, [[3,3,3]])
        kpoints = Kpoints.monkhorst_automatic([2,2,2], [0,0,0])
        self.assertEqual(kpoints.style, "Monkhorst")
        self.assertEqual(kpoints.kpts, [[2,2,2]])
        kpoints = Kpoints.automatic(100)
        self.assertEqual(kpoints.style, "Automatic")
        self.assertEqual(kpoints.kpts, [[100]])
        
class  PotcarTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','POTCAR')
        potcar = Potcar.from_file(filepath)
        self.assertEqual(potcar.symbols,["Fe","P","O"],"Wrong symbols read in for POTCAR")
        
    def test_potcar_map(self):
        fe_potcar = open(os.path.join(module_dir, 'vasp_testfiles','Fe_POTCAR')).read()
        #specify V instead of Fe - this makes sure the test won't pass if the code just grabs the POTCAR from the config file (the config file would grab the V POTCAR)
        potcar = Potcar(["V"], sym_potcar_map={"V": fe_potcar})
        self.assertEqual(potcar.symbols,["Fe"],"Wrong symbols read in for POTCAR")
        
class  VasprunTest(unittest.TestCase):
    
    def setUp(self):
        self.filepath = os.path.join(module_dir, 'vasp_testfiles','vasprun.xml')
        self.vasprun = Vasprun(self.filepath)
    
    def test_properties(self):
        
        vasprun = self.vasprun
        
        
        totalscsteps = sum([len(i['electronic_steps']) for i in vasprun.ionic_steps])
                
        self.assertEquals(29, len(vasprun.ionic_steps), "Incorrect number of energies read from vasprun.xml")
        self.assertEquals(308, totalscsteps, "Incorrect number of energies read from vasprun.xml")
        
        
        self.assertEquals([u'Li', u'Fe', u'Fe', u'Fe', u'Fe', u'P', u'P', u'P', u'P', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O', u'O']
, vasprun.atomic_symbols, "Incorrect symbols read from vasprun.xml")
        self.assertEquals(vasprun.final_structure.composition.reduced_formula, "LiFe4(PO4)4", "Wrong formula for final structure read.")
        self.assertIsNotNone(vasprun.incar, "Incar cannot be read")
        self.assertIsNotNone(vasprun.kpoints, "Kpoints cannot be read")
        self.assertIsNotNone(vasprun.eigenvalues, "Eigenvalues cannot be read")
        self.assertAlmostEqual(vasprun.final_energy,-269.38319884,7, "Wrong final energy")
        self.assertAlmostEqual(vasprun.tdos.get_gap(),2.0589,4,"Wrong gap from dos!")
        
        expectedans = (2.539, 4.0906, 1.5516, False)
        (gap, cbm, vbm, direct) = vasprun.eigenvalue_band_properties
        self.assertAlmostEqual(gap, expectedans[0]) 
        self.assertAlmostEqual(cbm, expectedans[1])
        self.assertAlmostEqual(vbm, expectedans[2])
        self.assertEqual(direct, expectedans[3])
        
        self.assertEqual(vasprun.potcar_symbols, [u'PAW_PBE Li 17Jan2003', u'PAW_PBE Fe 06Sep2000', u'PAW_PBE Fe 06Sep2000', u'PAW_PBE P 17Jan2003', u'PAW_PBE O 08Apr2002'])
        self.assertIsNotNone(vasprun.kpoints, "Kpoints cannot be read")
        self.assertIsNotNone(vasprun.actual_kpoints, "Actual kpoints cannot be read")
        self.assertIsNotNone(vasprun.actual_kpoints_weights, "Actual kpoints weights cannot be read")
        for atomdoses in vasprun.pdos:
            for orbitaldos in atomdoses:
                self.assertIsNotNone(orbitaldos, "Partial Dos cannot be read")
            
class OutcarTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','OUTCAR')
        outcar = Outcar(filepath)
        expected_mag = ({'d': 0.0, 'p': 0.003, 's': 0.002, 'tot': 0.005},
 {'d': 0.798, 'p': 0.008, 's': 0.007, 'tot': 0.813},
 {'d': 0.798, 'p': 0.008, 's': 0.007, 'tot': 0.813},
 {'d': 0.0, 'p': -0.117, 's': 0.005, 'tot': -0.112},
 {'d': 0.0, 'p': -0.165, 's': 0.004, 'tot': -0.162},
 {'d': 0.0, 'p': -0.117, 's': 0.005, 'tot': -0.112},
 {'d': 0.0, 'p': -0.165, 's': 0.004, 'tot': -0.162})
        expected_chg = ({'p': 0.154, 's': 0.078, 'd': 0.0, 'tot': 0.232}, {'p': 0.707, 's': 0.463, 'd': 8.316, 'tot': 9.486}, {'p': 0.707, 's': 0.463, 'd': 8.316, 'tot': 9.486}, {'p': 3.388, 's': 1.576, 'd': 0.0, 'tot': 4.964}, {'p': 3.365, 's': 1.582, 'd': 0.0, 'tot': 4.947}, {'p': 3.388, 's': 1.576, 'd': 0.0, 'tot': 4.964}, {'p': 3.365, 's': 1.582, 'd': 0.0, 'tot': 4.947})

        self.assertAlmostEqual(outcar.magnetization, expected_mag, 5, "Wrong magnetization read from Outcar")
        self.assertAlmostEqual(outcar.charge, expected_chg, 5, "Wrong charge read from Outcar")
        self.assertFalse(outcar.is_stopped)
        self.assertEqual(outcar.run_stats, {'System time (sec)': 0.938, 'Total CPU time used (sec)': 545.142, 'Elapsed time (sec)': 546.709, 'Maximum memory used (kb)': 0.0, 'Average memory used (kb)': 0.0, 'User time (sec)': 544.204})
        filepath = os.path.join(module_dir, 'vasp_testfiles','OUTCAR.stopped')
        outcar = Outcar(filepath)
        self.assertTrue(outcar.is_stopped)

class OszicarTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','OSZICAR')
        oszicar = Oszicar(filepath)
        self.assertEqual(len(oszicar.electronic_steps), len(oszicar.ionic_steps))
        self.assertEqual(len(oszicar.all_energies), 60)
        self.assertAlmostEqual(oszicar.final_energy, -526.63928)
        
if __name__ == '__main__':
    unittest.main()

