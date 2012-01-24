#!/usr/bin/env python
from __future__ import division

'''
Created on Sep 23, 2011
'''

__author__="Shyue Ping Ong"
__copyright__ = "Copyright 2011, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyue@mit.edu"
__date__ = "Sep 23, 2011"

import unittest

from pymatgen.transformations.standard_transformations import *
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.ewald import EwaldSummation

class TransformationsTest(unittest.TestCase):
    
    def setUp(self):
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.5,0.75])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        self.struct = Structure(lattice,["Si"] *2 ,coords)

    def test_identity_transformation(self):
        t = IdentityTransformation()
        self.assertEqual(self.struct, t.apply_transformation(self.struct))
    
    def test_rotation_transformation(self):
        t = RotationTransformation([0,1,0], 30, False)
        s2 = t.apply_transformation(self.struct)
        s1 = t.inverse.apply_transformation(s2)
        self.assertTrue((abs(s1.lattice.matrix - self.struct.lattice.matrix) < 1e-8).all())


class RemoveSpeciesTransformationTest(unittest.TestCase):
    
    def test_apply_transformation(self):
        t = RemoveSpeciesTransformation(["Li+"])
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.75,0.75])
        coords.append([0.5,0.5,0.5])
        coords.append([0.25,0.25,0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        struct = Structure(lattice,["Li+", "Li+", "O2-", "O2-"],coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s.composition.formula, "O2")
        

class SubstitutionTransformationTest(unittest.TestCase):
    
    def test_apply_transformation(self):
        t = SubstitutionTransformation({"Li+":"Na+", "O2-":"S2-"})
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.75,0.75])
        coords.append([0.5,0.5,0.5])
        coords.append([0.25,0.25,0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        struct = Structure(lattice,["Li+", "Li+", "O2-", "O2-"],coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s.composition.formula, "Na2 S2")

class SupercellTransformationTest(unittest.TestCase):
    
    def test_apply_transformation(self):
        t = SupercellTransformation([[2,1,0],[0,2,0],[1,0,2]])
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.75,0.75])
        coords.append([0.5,0.5,0.5])
        coords.append([0.25,0.25,0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        struct = Structure(lattice,["Li+", "Li+", "O2-", "O2-"],coords)
        s = t.apply_transformation(struct)
        self.assertEqual(s.composition.formula, "Li16 O16")

class PartialRemoveSpecieTransformationTest(unittest.TestCase):
    
    def test_apply_transformation(self):
        t = PartialRemoveSpecieTransformation("Li+", 1.0/3)
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.75,0.75])
        coords.append([0.5,0.5,0.5])
        coords.append([0.25,0.25,0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        struct = Structure(lattice,["Li+", "Li+", "Li+", "O2-"],coords)
        t.apply_transformation(struct)
        self.assertEqual(len(t.all_structures), 3)
        
class PartialRemoveSpecieTransformation2Test(unittest.TestCase):
    
    def test_apply_transformation(self):
        t = PartialRemoveSpecieTransformation2("Li+", 1.0/3)
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.75,0.75])
        coords.append([0.5,0.5,0.5])
        coords.append([0.25,0.25,0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        struct = Structure(lattice,["Li+", "Li+", "Li+", "O2-"],coords)
        self.assertAlmostEqual(EwaldSummation(t.apply_transformation(struct)).total_energy, -32.3858881789, 4, 'wrong structure chosen')
        

class OrderDisorderedStructureTransformationTest(unittest.TestCase):

    def test_apply_transformation(self):
        t = OrderDisorderedStructureTransformation()
        coords = list()
        coords.append([0,0,0])
        coords.append([0.75,0.75,0.75])
        coords.append([0.5,0.5,0.5])
        coords.append([0.25,0.25,0.25])
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        struct = Structure(lattice,[{"Si4+":0.5, "O2-": 0.25, "P5+": 0.25}, {"Si4+":0.5, "O2-": 0.25, "P5+": 0.25}, {"Si4+":0.5, "O2-": 0.25, "P5+": 0.25}, {"Si4+":0.5, "O2-": 0.25, "P5+": 0.25}] ,coords)
        t.apply_transformation(struct)
        self.assertEqual(len(t.all_structures), 12)
        
        struct = Structure(lattice,[{"Si4+":0.5}, {"Si4+":0.5}, {"P5+":0.5, "O2-": 0.5}, {"P5+":0.5, "O2-": 0.5}] ,coords)
        t.apply_transformation(struct)
        self.assertEqual(len(t.all_structures), 4)
        
        struct = Structure(lattice,[{"Si4+":0.333}, {"Si4+":0.333}, {"Si4+":0.333}, "O2-"] ,coords)
        t.apply_transformation(struct)
        self.assertEqual(len(t.all_structures), 3)
        
class PrimitiveCellTransformationTest(unittest.TestCase):
    
    def test_apply_transformation(self):
        t = PrimitiveCellTransformation()
        coords = list()
        coords.append([0,0,0])
        coords.append([0.375,0.375,0.375])
        coords.append([.5,.5,.5])
        coords.append([0.875,0.875,0.875])
        coords.append([0.125,0.125,0.125])
        coords.append([0.25,0.25,0.25])
        coords.append([0.625,0.625,0.625])
        coords.append([0.75,0.75,0.75])
        
        lattice = Lattice([[ 3.8401979337, 0.00, 0.00],[1.9200989668, 3.3257101909, 0.00],[0.00,-2.2171384943,3.1355090603]])
        struct = Structure(lattice,["Li+", "Li+","Li+", "Li+", "O2-", "O2-", "O2-", "O2-"],coords)
        s = t.apply_transformation(struct)
        self.assertEqual(len(s), 4)
        
class TransformationJsonTest(unittest.TestCase):
    
    def test_from_json(self):
        self.assertIsInstance(transformation_from_json('{"name": "IdentityTransformation", "init_args": {}}'), IdentityTransformation)
        self.assertIsInstance(transformation_from_json('{"name": "RotationTransformation", "init_args": {"angle": 30, "angle_in_radians": false, "axis": [0, 1, 0]}}'), RotationTransformation)
        
if __name__ == "__main__":
    unittest.main()