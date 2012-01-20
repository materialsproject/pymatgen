#!/usr/bin/python
import unittest
import os

from pymatgen.io.vaspio_set import MITVaspInputSet, MaterialsProjectVaspInputSet
from pymatgen.io.vaspio import Poscar
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from numpy import array

module_dir = os.path.dirname(os.path.abspath(__file__))

class MITMaterialsProjectVaspInputSetTest(unittest.TestCase):
    
    def setUp(self):
        self.mitparamset = MITVaspInputSet()
        filepath = os.path.join(module_dir,'vasp_testfiles','POSCAR')
        poscar = Poscar.from_file(filepath)
        self.struct = poscar.struct
        
        self.paramset = MaterialsProjectVaspInputSet()
        filepath = os.path.join(module_dir,'vasp_testfiles','POSCAR')
        poscar = Poscar.from_file(filepath)
        self.struct = poscar.struct
    
    def test_get_potcar_symbols(self):
        syms = self.paramset.get_potcar_symbols(self.struct)
        self.assertEquals(syms, ['Fe_pv', 'P', 'O'])
        
    def test_get_incar(self):
        incar = self.paramset.get_incar(self.struct)
        self.assertEqual(incar['LDAUU'], [5.3, 0, 0])
        self.assertAlmostEqual(incar['EDIFF'], 0.0012)
        
        incar = self.mitparamset.get_incar(self.struct)
        self.assertEqual(incar['LDAUU'], [4.0, 0, 0])
        self.assertAlmostEqual(incar['EDIFF'], 0.0012)
         
        si = 14
        coords = list()
        coords.append(array([0,0,0]))
        coords.append(array([0.75,0.5,0.75]))

        #Silicon structure for testing.
        latt = Lattice(array([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]]))
        struct = Structure(latt,[si,si],coords)
        incar = incar = self.paramset.get_incar(struct)
        self.assertNotIn("LDAU", incar)
        
    def test_get_kpoints(self):
        kpoints = self.paramset.get_kpoints(self.struct)
        self.assertEquals(kpoints.kpts, [[2,4,6]])
        self.assertEquals(kpoints.style, 'Monk')
        
        kpoints = self.mitparamset.get_kpoints(self.struct)
        self.assertEquals(kpoints.kpts, [[2,4,4]])
        self.assertEquals(kpoints.style, 'Monk')
        
if __name__ == '__main__':
    unittest.main()

