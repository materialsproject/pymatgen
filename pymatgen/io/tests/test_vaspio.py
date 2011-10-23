#!/usr/bin/python
import unittest
import os

from pymatgen.io.vaspio import Poscar, Potcar, Kpoints, Incar, Vasprun, Outcar
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Composition, Structure
from numpy import array

module_dir = os.path.dirname(os.path.abspath(__file__))

class  PoscarTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','POSCAR.gz')
        poscar = Poscar.from_file(filepath)
        comp = poscar.struct.composition
        self.assertEqual(comp,Composition.from_formula("Fe4P4O16"))

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
        filepath = os.path.join(module_dir, 'vasp_testfiles','INCAR.gz')
        incar = Incar.from_file(filepath)
        incar["LDAU"] = "T"
        self.assertEqual(incar["ALGO"],"Damped","Wrong Algo")
        self.assertEqual(float(incar["EDIFF"]),1e-4,"Wrong EDIFF")
        
    def test_diff(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','INCAR.gz')
        incar1 = Incar.from_file(filepath)
        filepath = os.path.join(module_dir, 'vasp_testfiles','INCAR.2.gz')
        incar2 = Incar.from_file(filepath)
                
        self.assertEqual(incar1.diff(incar2), {'Same parameters': {'IBRION': 2, 'PREC': 'Accurate', 'ISIF': 3, 'LMAXMIX': 4, 'LREAL': 'Auto', 'ISPIN': 2, 'LORBIT': '11', 'SIGMA': 0.05}, 'Different': {'MAGMOM': {'INCAR1': [6, -6, -6, 6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6], 'INCAR2': 'Default'}, 'NKRED': {'INCAR1': '2', 'INCAR2': 'Default'}, 'ENCUTFOCK': {'INCAR1': '0', 'INCAR2': 'Default'}, 'NUPDOWN': {'INCAR1': '0', 'INCAR2': 'Default'}, 'EDIFF': {'INCAR1': '1E-4', 'INCAR2': 0.0001}, 'HFSCREEN': {'INCAR1': '0.207', 'INCAR2': 'Default'}, 'LSCALU': {'INCAR1': '.FALSE.', 'INCAR2': 'Default'}, 'SYSTEM': {'INCAR1': 'id=[0] dblock_code=[97763-ICSD] formula=[Li Mn (P O4)] sg_name=[P n m a]', 'INCAR2': 'id=[91090] dblock_code=[20070929235612LiNiO-59.53134651-VASP] formula=[Li3 Ni3 O6] sg_name=[R-3m]'}, 'ENCUT': {'INCAR1': '500', 'INCAR2': 'Default'}, 'NSIM': {'INCAR1': '1', 'INCAR2': 'Default'}, 'LCHARG': {'INCAR1': '.TRUE.', 'INCAR2': 'Default'}, 'LPLANE': {'INCAR1': '.TRUE.', 'INCAR2': 'Default'}, 'ALGO': {'INCAR1': 'Damped', 'INCAR2': 'Fast'}, 'LHFCALC': {'INCAR1': '.TRUE.', 'INCAR2': 'Default'}, 'TIME': {'INCAR1': '0.4', 'INCAR2': 'Default'}, 'ISMEAR': {'INCAR1': 0, 'INCAR2': -5}, 'LWAVE': {'INCAR1': True, 'INCAR2': False}, 'NPAR': {'INCAR1': 8, 'INCAR2': 1}, 'NSW': {'INCAR1': 99, 'INCAR2': 51}, 'ISPIND': {'INCAR1': '2', 'INCAR2': 'Default'}}})

class  KpointsTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','KPOINTS.auto')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts,[[10]],"Wrong kpoint lattice read")
        filepath = os.path.join(module_dir, 'vasp_testfiles','KPOINTS.gz')
        kpoints = Kpoints.from_file(filepath)
        self.assertEqual(kpoints.kpts,[[2,4,6]],"Wrong kpoint lattice read")

class  PotcarTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','POTCAR.gz')
        potcar = Potcar.from_file(filepath)
        self.assertEqual(potcar.symbols,["Fe","P","O"],"Wrong symbols read in for POTCAR")
        
class  VasprunTest(unittest.TestCase):
    
    def setUp(self):
        self.filepath = os.path.join(module_dir, 'vasp_testfiles','vasprun.xml.gz')
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
        print vasprun.kpoints
        print vasprun.kpoints.to_dict
        self.assertIsNotNone(vasprun.eigenvalues, "Eigenvalues cannot be read")
        self.assertAlmostEqual(vasprun.final_energy,-269.38319884,7, "Wrong final energy")
        self.assertAlmostEqual(vasprun.tdos.get_gap(),2.0589,4,"Wrong gap from dos!")
        self.assertEqual(vasprun.potcar_symbols, [u'PAW_PBE Li 17Jan2003', u'PAW_PBE Fe 06Sep2000', u'PAW_PBE Fe 06Sep2000', u'PAW_PBE P 17Jan2003', u'PAW_PBE O 08Apr2002'])
        self.assertIsNotNone(vasprun.kpoints, "Kpoints cannot be read")
        self.assertIsNotNone(vasprun.actual_kpoints, "Actual kpoints cannot be read")
        self.assertIsNotNone(vasprun.actual_kpoints_weights, "Actual kpoints weights cannot be read")
        for atomdoses in vasprun.pdos:
            for orbitaldos in atomdoses:
                self.assertIsNotNone(orbitaldos, "Partial Dos cannot be read")
            
class  OutcarTest(unittest.TestCase):
    
    def test_init(self):
        filepath = os.path.join(module_dir, 'vasp_testfiles','OUTCAR.gz')
        outcar = Outcar(filepath)
        expected_ans = {1: 0.005, 2: 0.813, 3: 0.813, 4: -0.112, 5: -0.162, 6: -0.112, 7: -0.162}
        self.assertAlmostEqual(outcar.get_magnetization(), expected_ans, 5, "Wrong magnetization read from Outcar")

if __name__ == '__main__':
    unittest.main()

