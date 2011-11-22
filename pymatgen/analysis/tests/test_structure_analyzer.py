import unittest
import os

from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder, solid_angle
import pymatgen.io.vaspio 
from pymatgen.core.periodic_table import Element

class VoronoiCoordFinderTest(unittest.TestCase):

    def setUp(self):
        module_path = os.path.dirname(pymatgen.io.vaspio.__file__)
        filepath = os.path.join(module_path, 'tests','vasp_testfiles', 'vasprun.xml')
        reader = pymatgen.io.vaspio.Vasprun(filepath)
        s = reader.final_structure
        self.finder = VoronoiCoordFinder(s,[Element("O")])
    
    def test_get_voronoi_polyhedra(self):
        self.assertEqual(len(self.finder.get_voronoi_polyhedra(0).items()),10, "Incorrect number of results returned for get_voronoi_polyhedra")
        
    def test_get_coordination_number(self):
        print self.finder.get_coordination_number(0)
        self.assertAlmostEqual(self.finder.get_coordination_number(0), 5.60588600732, 7, "Incorrect coordination number returned!")                

    def test_get_coordinated_sites(self):
        self.assertEqual(len(self.finder.get_coordinated_sites(0)), 10)

class MiscFunctionTest(unittest.TestCase):

    def test_solid_angle(self):
        center = [2.294508207929496, 4.4078057081404, 2.299997773791287] 
        coords = [[1.627286218099362,3.081185538926995,3.278749383217061],[1.776793751092763,2.93741167455471,3.058701096568852],[3.318412187495734,2.997331084033472,2.022167590167672],[3.874524708023352,4.425301459451914,2.771990305592935],[2.055778446743566,4.437449313863041,4.061046832034642]]
        self.assertAlmostEqual(solid_angle(center, coords), 1.83570965938, 7, "Wrong result returned by solid_angle")
    
if __name__ == '__main__':
    unittest.main()
