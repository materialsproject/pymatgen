import unittest

from pymatgen.core.structure import Structure
import subprocess
import os 
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io import atat
from pymatgen.util.testing import PymatgenTest
from pymatgen.command_line.mcsqs_caller import run_mcsqs
from monty.os.path import which

__author__ = "Handong ling"
__version__ = "0.1"
__maintainer__ = "Handong Ling"
__email__ = "handongling@berkeley.edu"
__date__ = "June 2019"

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        'test_files', "mcsqs")

@unittest.skipIf(not which('mcsqs'), "mcsqs executable not present")
class Mcsqs_CallerTest(PymatgenTest):
	def setUp(self):
		self.pztstrings = np.load(os.path.join(test_dir, "pztstrings.npy"), allow_pickle=True)
		self.pztstrings2 = np.load(os.path.join(test_dir, "pztstrings2.npy"), allow_pickle=True)
		self.struc = self.get_structure('Pb2TiZrO6')
		self.struc.replace_species({'Ti': {'Ti':0.5, 'Zr':0.5 }, 'Zr': {'Ti':0.5, 'Zr':0.5 }})


	def test_Mcsqs_Caller_supercell(self):
		count = 0
		strings = []
		passed = [] 

		sqs = run_mcsqs(self.struc, {2:6, 3:4}, supercell = [2,1,1], total_atoms = None, max_time = None)
		self.assertEqual(atat.Mcsqs(sqs).to_string() in self.pztstrings, True)
		os.remove('sqscell.out')
		os.remove('rndstrgrp.out')
		os.remove('bestcorr.out')
		os.remove('rndstr.in')
		os.remove('sym.out')
		os.remove('mcsqs.log')
		os.remove('bestsqs.out')
		os.remove('clusters.out')

	def test_Mcsqs_Caller_total_atoms(self):
		sqs = run_mcsqs(self.struc, {2:6, 3:4}, total_atoms = 20, max_time = None)
		self.assertEqual(atat.Mcsqs(sqs).to_string() in self.pztstrings2, True)
		os.remove('sqscell.out')
		os.remove('rndstrgrp.out')
		os.remove('bestcorr.out')
		os.remove('rndstr.in')
		os.remove('sym.out')
		os.remove('mcsqs.log')
		os.remove('bestsqs.out')
		os.remove('clusters.out')

	def test_Mcsqs_Caller_supercell_timeout(self):
		# sqs = run_mcsqs(self.struc, {2:6, 3:4}, supercell = 1000, max_time = 0.000001)
		self.assertRaises(TimeoutError, run_mcsqs, self.struc, {2:6, 3:4}, supercell = [2,1,1], max_time = 0.000001)		
		# os.remove('sqscell.out')
		# os.remove('rndstrgrp.out')
		# os.remove('bestcorr.out')
		# os.remove('rndstr.in')
		# os.remove('sym.out')
		# os.remove('mcsqs.log')
		# os.remove('clusters.out')
