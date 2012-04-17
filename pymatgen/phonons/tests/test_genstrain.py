import unittest
import os, sys
sys.path.append('/home/MDEJONG1/pythonplayground/pymatgen/pymatgen_repo/pymatgen_repo')
from pymatgen.io.cifio import CifParser


#print CifParser('aluminum.cif').get_structures(False)[0]

class GenstrainTest(unittest.TestCase):

	def setUp(self):
		self.structure = CifParser(os.path.join(self.module_dir, 'tests', 'aluminum.cif')).get_structures(False)[0]
		
	
	def test_deformation(self):




#		print self.structure
#		print 'fjwk'

if __name__=='__main__':
	unittest.main()

#A = GenstrainTest.setUp()
#print type(A)

